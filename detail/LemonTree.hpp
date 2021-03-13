/**
 * @file LemonTree.hpp
 * @brief Implementation of the class for learning module networks
 *        using the approach of Lemon Tree.
 * @author Ankit Srivastava <asrivast@gatech.edu>
 *
 * Copyright 2020 Georgia Institute of Technology
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
#ifndef DETAIL_LEMONTREE_HPP_
#define DETAIL_LEMONTREE_HPP_

#include "Ganesh.hpp"
#include "Module.hpp"
#include "ConsensusCluster.hpp"

#include "mxx/distribution.hpp"

#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/gzip.hpp>


template <typename Data, typename Var, typename Set>
LemonTree<Data, Var, Set>::LemonTree(
  const mxx::comm& comm,
  const Data& data
) : ModuleNetworkLearning<Data, Var, Set>(comm, data)
{
  TIMER_RESET(m_tWrite);
  TIMER_RESET(m_tGanesh);
  TIMER_RESET(m_tConsensus);
  TIMER_RESET(m_tModules);
  TIMER_RESET(m_tSync);
}

template <typename Data, typename Var, typename Set>
/**
 * @brief Default destructor.
 */
LemonTree<Data, Var, Set>::~LemonTree(
)
{
  if (this->m_comm.is_first()) {
    TIMER_ELAPSED_NONZERO("Time taken in writing the files: ", m_tWrite);
    TIMER_ELAPSED_NONZERO("Time taken in the GaneSH run: ", m_tGanesh);
    TIMER_ELAPSED_NONZERO("Time taken in consensus clustering: ", m_tConsensus);
    TIMER_ELAPSED_NONZERO("Time taken in learning the modules: ", m_tModules);
    TIMER_ELAPSED_NONZERO("Time taken in synchronizing the modules: ", m_tSync);
  }
}

template <typename Data, typename Var, typename Set>
template <typename Generator>
std::list<Set>
LemonTree<Data, Var, Set>::singleGaneshRun(
  Generator& generator,
  const pt::ptree& ganeshConfigs
) const
{
  auto numSteps = ganeshConfigs.get<uint32_t>("num_steps");
  auto initClusters = ganeshConfigs.get<Var>("init_num_clust");
  if ((initClusters == 0) || (initClusters > this->m_data.numVars())) {
    initClusters = this->m_data.numVars() / 2;
  }
  Ganesh<Data, Var, Set> ganesh(this->m_data);
  ganesh.initializeRandom(generator, initClusters);
  for (auto s = 0u; s <= numSteps; ++s) {
    LOG_MESSAGE(info, "Step %u", s);
    ganesh.clusterTwoWay(generator, this->m_comm);
  }
  // XXX: Lemon Tree writes out only the last sampled cluster
  // and uses that for the downstream tasks per run
  // We can replicate this behavior by storing only the last sample in each run
  // Further, we do not use -burn_in and -sample_steps because they do not
  // have any effect on the outcome
  LOG_MESSAGE(info, "Sampling");
  const auto& primaryClusters = ganesh.primaryClusters();
  std::list<Set> varClusters;
  for (const auto& cluster : primaryClusters) {
    varClusters.push_back(cluster.elements());
  }
  return varClusters;
}

template <typename Data, typename Var, typename Set>
template <typename Generator>
std::list<std::list<Set>>
LemonTree<Data, Var, Set>::clusterVarsGanesh(
  const pt::ptree& ganeshConfigs
) const
{
  auto randomSeed = ganeshConfigs.get<uint32_t>("seed");
  auto numRuns = ganeshConfigs.get<uint32_t>("num_runs");
  Generator generator(randomSeed);
  std::list<std::list<Set>> sampledClusters;
  for (auto r = 0u; r < numRuns; ++r) {
    LOG_MESSAGE(info, "Run %u", r);
    // XXX: Seeding in this way to compare the results with Lemon-Tree;
    //      Otherwise, we can carry over generator state across runs
    generator.seed(randomSeed + r);
    auto varClusters = this->singleGaneshRun(generator, ganeshConfigs);
    sampledClusters.push_back(varClusters);
  }
  return sampledClusters;
}

template <typename Data, typename Var, typename Set>
void
LemonTree<Data, Var, Set>::writeVarClusters(
  const std::string& clusterFile,
  const std::list<std::list<Set>>& varClusters
) const
{
  LOG_MESSAGE(info, "Writing variable clusters to %s", clusterFile);
  std::ofstream gf(clusterFile);
  auto g = 0u;
  for (auto git = varClusters.begin(); git != varClusters.end(); ++git, ++g) {
    std::string clustersFile(clusterFile + "." + std::to_string(g));
    std::ofstream cf(clustersFile);
    auto c = 0u;
    for (auto cit = git->begin(); cit != git->end(); ++cit, ++c) {
      for (const auto e : *cit) {
        cf << this->m_data.varName(e) << "\t" << c << std::endl;
      }
    }
    gf << clustersFile << std::endl;
  }
}

template <typename Data, typename Var, typename Set>
std::multimap<Var, Var>
LemonTree<Data, Var, Set>::clusterConsensus(
  const std::list<std::list<Set>>&& varClusters,
  const pt::ptree& consensusConfigs
) const
{
  auto minWeight = consensusConfigs.get<double>("min_weight");
  auto tolerance = consensusConfigs.get<double>("tolerance", 1e-5);
  auto maxSteps = consensusConfigs.get<uint32_t>("max_steps", 1000);
  auto minClustSize = consensusConfigs.get<uint32_t>("min_clust_size");
  auto minClustScore = consensusConfigs.get<double>("min_clust_score");
  auto result = consensusCluster(std::move(varClusters), this->m_data.numVars(), minWeight,
                                 tolerance, maxSteps, minClustSize, minClustScore);
  return result;
}

template <typename Data, typename Var, typename Set>
void
LemonTree<Data, Var, Set>::writeConsensusCluster(
  const std::string& consensusFile,
  const std::multimap<Var, Var>& vertexClusters
) const
{
  LOG_MESSAGE(info, "Writing consensus clusters to %s", consensusFile);
  std::ofstream out(consensusFile);
  // Write output in tab seperated format
  for (const auto& cx : vertexClusters) {
    out << this->m_data.varName(cx.second) << "\t" << static_cast<uint32_t>(cx.first) << std::endl;
  }
}

template <typename Data, typename Var, typename Set>
template <typename Generator>
std::list<std::list<Set>>
LemonTree<Data, Var, Set>::clusterObsGanesh(
  const uint32_t numRuns,
  const uint32_t numSteps,
  const uint32_t burnSteps,
  const uint32_t sampleSteps,
  Generator& generator,
  const Set& clusterVars
) const
{
  std::list<std::list<Set>> sampledClusters;
  // Initialize Gibbs sampler algorithm for this cluster
  Ganesh<Data, Var, Set> ganesh(this->m_data);
  ganesh.initializeGiven(generator, std::list<Set>(1, clusterVars));
  for (auto r = 0u; r < numRuns; ++r) {
    auto s = 0u;
    for ( ; s < burnSteps; ++s) {
      LOG_MESSAGE(info, "Step %u (burn in)", s);
      ganesh.clusterSecondary(generator);
    }
    for ( ; s < numSteps; ++s) {
      LOG_MESSAGE(info, "Step %u (sampling)", s);
      ganesh.clusterSecondary(generator);
      if ((numSteps - (s + 1)) % sampleSteps == 0) {
        LOG_MESSAGE(info, "Sampling");
        // There should be only one primary cluster; get a reference to it
        const auto& primaryCluster = ganesh.primaryClusters().front();
        const auto& secondaryClusters = primaryCluster.secondaryClusters();
        std::list<Set> obsClusters;
        for (const auto& cluster : secondaryClusters) {
          obsClusters.push_back(cluster.elements());
        }
        sampledClusters.push_back(obsClusters);
      }
    }
  }
  return sampledClusters;
}

template <typename Data, typename Var, typename Set>
void
LemonTree<Data, Var, Set>::readCandidateParents(
  const std::string& fileName,
  Set& candidateParents
) const
{
  LOG_MESSAGE(info, "Reading candidate parents from %s", fileName);
  std::ifstream regFile(boost::filesystem::canonical(fileName).string());
  std::string name;
  while (std::getline(regFile, name)) {
    auto v = this->m_data.varIndex(name);
    if (v < this->m_data.numVars()) {
      candidateParents.insert(v);
    }
  }
  LOG_MESSAGE(info, "Read %u candidate parents", candidateParents.size());
}

template <typename Data, typename Var, typename Set>
template <typename Generator>
std::list<Module<Data, Var, Set>>
LemonTree<Data, Var, Set>::constructModulesWithTrees(
  const std::multimap<Var, Var>&& coClusters,
  Generator& generator,
  const pt::ptree& modulesConfigs
) const
{
  auto numRuns = modulesConfigs.get<uint32_t>("num_runs");
  auto numSteps = modulesConfigs.get<uint32_t>("num_steps");
  auto burnSteps = modulesConfigs.get<uint32_t>("burn_in");
  auto sampleSteps = modulesConfigs.get<uint32_t>("sample_steps");
  auto scoreBHC = modulesConfigs.get<bool>("use_bayesian_score", true);
  auto scoreGain = modulesConfigs.get<double>("score_gain");
  std::list<Module<Data, Var, Set>> modules;
  auto m = 0u;
  for (auto cit = coClusters.begin(); cit != coClusters.end(); ++m) {
    LOG_MESSAGE(info, "Module %u: Learning tree structures", m);
    // Get the range of variables in this cluster
    auto clusterIts = coClusters.equal_range(cit->first);
    // Add all the variables in the cluster to a set
    auto clusterVars = Set(this->m_data.numVars());
    for (auto vit = clusterIts.first; vit != clusterIts.second; ++vit) {
      clusterVars.insert(vit->second);
    }
    // Create a module for this variable cluster
    modules.emplace_back(std::move(clusterVars), this->m_comm, this->m_data);
    auto& module = modules.back();
    // Sample observation clusters for this module
    auto sampledClusters = this->clusterObsGanesh(numRuns, numSteps, burnSteps, sampleSteps,
                                                  generator, module.variables());
    // Learn tree structures from the observation clusters
    module.learnTreeStructures(std::move(sampledClusters), scoreBHC, scoreGain);
    cit = clusterIts.second;
  }
  return modules;
}

template <typename Data, typename Var, typename Set>
template <typename Generator>
void
LemonTree<Data, Var, Set>::learnModulesParents(
  std::list<Module<Data, Var, Set>>& modules,
  Generator& generator,
  const pt::ptree& modulesConfigs
) const
{
  auto regFile = modulesConfigs.get<std::string>("reg_file");
  auto betaMax = modulesConfigs.get<double>("beta_reg");
  auto numSplits = modulesConfigs.get<uint32_t>("num_reg");
  Set candidateParents(this->m_data.numVars());
  if (!regFile.empty()) {
    // Read candidate parents from the given file
    this->readCandidateParents(regFile, candidateParents);
  }
  else {
    // Add all the variables as candidate parents
    for (Var v = 0u; v < candidateParents.max(); ++v) {
      candidateParents.insert(v);
    }
  }
  OptimalBeta ob(0.0, betaMax, 1e-5);
  auto m = 0u;
  for (auto moduleIt = modules.begin(); moduleIt != modules.end(); ++moduleIt, ++m) {
    LOG_MESSAGE(info, "Module %u: Learning parents", m);
    moduleIt->learnParents(generator, candidateParents, ob, numSplits);
  }
  LOG_MESSAGE(info, "Done learning module parents");
}

template <typename Data, typename Var, typename Set>
template <typename Generator>
void
LemonTree<Data, Var, Set>::learnModulesParents_nodes(
  std::list<Module<Data, Var, Set>>& modules,
  Generator& generator,
  const Set&& candidateParents,
  const double betaMax,
  const uint32_t numSplits
) const
{
  OptimalBeta ob(0.0, betaMax, 1e-5);
  std::vector<uint32_t> moduleNodeCount(modules.size());
  std::vector<uint32_t> moduleNodeWeights;
  auto totalNodes = 0u;
  auto moduleCit = modules.cbegin();
  for (auto m = 0u; m < modules.size(); ++m, ++moduleCit) {
    moduleNodeCount[m] = moduleCit->nodeCount();
    auto nodeWeights = moduleCit->nodeWeights();
    moduleNodeWeights.insert(moduleNodeWeights.end(), nodeWeights.begin(), nodeWeights.end());
    totalNodes += moduleNodeCount[m];
  }
  std::vector<uint32_t> moduleNodeCountPrefix(moduleNodeCount.size());
  std::partial_sum(moduleNodeCount.cbegin(), moduleNodeCount.cend(), moduleNodeCountPrefix.begin());
  std::vector<uint32_t> moduleNodeWeightsPrefix(moduleNodeWeights.size());
  std::partial_sum(moduleNodeWeights.cbegin(), moduleNodeWeights.cend(), moduleNodeWeightsPrefix.begin());
  auto totalWeight = moduleNodeWeightsPrefix.back();
  mxx::blk_dist block(totalWeight, this->m_comm.size(), this->m_comm.rank());
  auto myFirstNode = std::distance(moduleNodeWeightsPrefix.cbegin(),
                                   std::lower_bound(moduleNodeWeightsPrefix.cbegin(),
                                                    moduleNodeWeightsPrefix.cend(),
                                                    block.eprefix_size()));
  if (!this->m_comm.is_first()) {
    myFirstNode += 1;
  }
  auto myLastNode = std::distance(moduleNodeWeightsPrefix.cbegin(),
                                  std::lower_bound(moduleNodeWeightsPrefix.cbegin(),
                                                   moduleNodeWeightsPrefix.cend(),
                                                   block.iprefix_size())) + 1;
  auto myNodeCount = myLastNode - myFirstNode;
  // First, advance the PRNG state to account for the number of
  // random generations for nodes on the previous ranks
  advance(generator, myFirstNode * 2 * numSplits);
  // XXX: We need to track the validity of nodes because some nodes may not learn any splits
  //      and that information will be required for synchornization later
  std::vector<uint8_t> myValidNodes(myNodeCount, 0);
  std::vector<std::tuple<Var, Var, double>> myNodeSplits(myNodeCount * 2 * numSplits);
  // Find the indices for the modules which contain the first and last node on this rank
  auto myFirstModule = std::distance(moduleNodeCountPrefix.cbegin(),
                                     std::lower_bound(moduleNodeCountPrefix.cbegin(),
                                                      moduleNodeCountPrefix.cend(),
                                                      myFirstNode + 1));
  auto myLastModule = std::distance(moduleNodeCountPrefix.cbegin(),
                                    std::lower_bound(moduleNodeCountPrefix.cbegin(),
                                                     moduleNodeCountPrefix.cend(),
                                                     myLastNode));
  auto moduleIt = std::next(modules.begin(), myFirstModule);
  auto validIt = myValidNodes.begin();
  auto splitIt = myNodeSplits.begin();
  auto prevNodeCount = myFirstNode;
  for (auto m = myFirstModule; m <= myLastModule; ++m, ++moduleIt) {
    auto firstNode = prevNodeCount - (moduleNodeCountPrefix[m] - moduleNodeCount[m]);
    auto numNodes = std::min(static_cast<uint32_t>(myLastNode), moduleNodeCountPrefix[m]) - prevNodeCount;
    if (numNodes > 0) {
      LOG_MESSAGE(info, "Module %u: Learning parents for %u nodes (starting node index = %u)", m, numNodes, firstNode);
      moduleIt->learnParents(generator, candidateParents, ob, numSplits, firstNode, numNodes, validIt, splitIt);
      prevNodeCount += numNodes;
    }
  }
  myNodeSplits.resize(std::distance(myNodeSplits.begin(), splitIt));
  LOG_MESSAGE(info, "Done learning module parents");
  this->m_comm.barrier();
  TIMER_START(m_tSync);
  // XXX: We can determine the number of elements on each processor for these calls
  //      and eliminate a gather call. Need to evaluate if it is worth the effort
  auto allValidNodes = mxx::allgatherv(myValidNodes, this->m_comm);
  auto allNodeSplits = mxx::allgatherv(myNodeSplits, this->m_comm);
  moduleIt = modules.begin();
  auto validCit = allValidNodes.cbegin();
  auto splitCit = allNodeSplits.cbegin();
  for (auto m = 0u; m < modules.size(); ++m, ++moduleIt) {
    // Synchronize the node parents for this module
    LOG_MESSAGE(info, "Module %u: Synchronizing parents for all nodes", m);
    moduleIt->syncParents(numSplits, validCit, splitCit);
  }
  LOG_MESSAGE(info, "Done synchronizing module parents");
  this->m_comm.barrier();
  TIMER_PAUSE(m_tSync);
}

template <typename Data, typename Var, typename Set>
template <typename Generator>
void
LemonTree<Data, Var, Set>::learnModulesParents_splits(
  std::list<Module<Data, Var, Set>>& modules,
  Generator& generator,
  const Set&& candidateParents,
  const double betaMax,
  const uint32_t numSplits
) const
{
  TIMER_DECLARE(tCandidates);
  OptimalBeta ob(0.0, betaMax, 1e-5);
  std::vector<uint32_t> moduleNodeCount(modules.size());
  std::vector<uint64_t> moduleMaxSplits(modules.size());
  auto moduleCit = modules.cbegin();
  for (auto m = 0u; m < modules.size(); ++m, ++moduleCit) {
    moduleNodeCount[m] = moduleCit->nodeCount();
    moduleMaxSplits[m] = moduleCit->maxSplits(candidateParents);
  }
  std::vector<uint64_t> moduleNodeCountPrefix(moduleNodeCount.size());
  std::partial_sum(moduleNodeCount.cbegin(), moduleNodeCount.cend(), moduleNodeCountPrefix.begin());
  std::vector<uint64_t> moduleMaxSplitsPrefix(moduleMaxSplits.size());
  std::partial_sum(moduleMaxSplits.cbegin(), moduleMaxSplits.cend(), moduleMaxSplitsPrefix.begin());
  auto totalSplits = moduleMaxSplitsPrefix.back();
  mxx::blk_dist block(totalSplits, this->m_comm.size(), this->m_comm.rank());
  auto myFirstModule = std::distance(moduleMaxSplitsPrefix.cbegin(),
                                     std::lower_bound(moduleMaxSplitsPrefix.cbegin(),
                                                      moduleMaxSplitsPrefix.cend(),
                                                      block.eprefix_size() + 1));
  auto myLastModule = std::distance(moduleMaxSplitsPrefix.cbegin(),
                                    std::lower_bound(moduleMaxSplitsPrefix.cbegin(),
                                                     moduleMaxSplitsPrefix.cend(),
                                                     block.iprefix_size()));
  std::vector<std::tuple<uint32_t, Var, Var, double>> mySplits;
  auto moduleIt = std::next(modules.begin(), myFirstModule);
  auto prevSplits = block.eprefix_size();
  for (auto m = myFirstModule; m <= myLastModule; ++m, ++moduleIt) {
    LOG_MESSAGE(info, "Module %u: Learning candidate parents splits", m);
    auto firstNode = moduleNodeCountPrefix[m] - moduleNodeCount[m];
    auto firstSplit = prevSplits - (moduleMaxSplitsPrefix[m] - moduleMaxSplits[m]);
    auto maxSplits = std::min(block.iprefix_size(), moduleMaxSplitsPrefix[m]) - prevSplits;
    if (maxSplits > 0) {
      moduleIt->candidateParentsSplits(mySplits, candidateParents, ob, firstNode, firstSplit, maxSplits);
      prevSplits += maxSplits;
    }
  }
  // Redistribute to get the same number of candidate splits on every processor
  mxx::stable_distribute_inplace(mySplits, this->m_comm);
  this->m_comm.barrier();
  if (this->m_comm.is_first()) {
    TIMER_ELAPSED("Time taken in learning candidate splits: ", tCandidates);
  }
  TIMER_DECLARE(tChoose);
  // Compute the max score for each node across all the processors
  std::vector<double> mySplitsScoresMax(moduleNodeCountPrefix.back(), std::numeric_limits<double>::lowest());
  std::for_each(mySplits.cbegin(), mySplits.cend(),
                [&mySplitsScoresMax] (const std::tuple<uint32_t, Var, Var, double>& split)
                                     { mySplitsScoresMax[std::get<0>(split)] = std::max(mySplitsScoresMax[std::get<0>(split)],
                                                                                        std::get<3>(split)); });
  auto allSplitsScoresMax = mxx::allreduce(mySplitsScoresMax, mxx::max<double>(), this->m_comm);
  // Set of local node indices
  // XXX: Using set instead of unordered_set because we want sorted indices
  std::set<uint32_t> myNodeIdx;
  // Count of node splits on this processor for all the nodes
  std::vector<uint64_t> mySplitsCounts(moduleNodeCountPrefix.back(), 0);
  // Vector of <node, weight> for each split, where weight = exp(score)
  std::vector<std::pair<uint32_t, double>> mySplitsWeights(mySplits.size());
  // Sum of weights of splits on this processor for all the nodes
  std::vector<double> mySplitsWeightsSum(moduleNodeCountPrefix.back(), 0);
  auto s = 0u;
  auto splitFirst = mySplits.cbegin();
  while (splitFirst != mySplits.cend()) {
    auto n = std::get<0>(*splitFirst);
    myNodeIdx.insert(n);
    auto splitLast = std::find_if(splitFirst, mySplits.cend(),
                                  [&n] (const std::tuple<uint32_t, Var, Var, double>& split)
                                       { return std::get<0>(split) != n; });
    mySplitsCounts[n] += std::distance(splitFirst, splitLast);
    auto nodeMaxScore = allSplitsScoresMax[n];
    for (auto it = splitFirst; it != splitLast; ++it, ++s) {
      auto weight = exp(std::get<3>(*it) - nodeMaxScore);
      mySplitsWeights[s] = std::make_pair(n, weight);
      mySplitsWeightsSum[n] += weight;
    }
    splitFirst = splitLast;
  }
  // Get the count of splits of the first node on previous processor(s)
  auto myLastNode = std::get<0>(mySplits.back());
  auto lastSplitsCount = std::make_pair(myLastNode, mySplitsCounts[myLastNode]);
  auto addSplitsCount = [] (const std::pair<uint32_t, uint64_t>& a,
                            const std::pair<uint32_t, uint64_t>& b)
                           { return (a.first == b.first) ?
                                    std::make_pair(b.first, a.second + b.second) : b; };
  auto firstCountsPrefix = mxx::exscan(lastSplitsCount, addSplitsCount, this->m_comm, false);
  // Get the total count of splits for all the nodes
  // Also get the sum of the weights of splits for each node on all the processors
  // XXX: We really only need these for the nodes on this processor
  //      Can we make it more efficient?
  auto allSplitsCounts = mxx::allreduce(mySplitsCounts, this->m_comm);
  auto allSplitsWeightsSum = mxx::allreduce(mySplitsWeightsSum, this->m_comm);
  // First, handle any nodes with infinite weights
  static auto findInf = [] (const double w) { return std::isinf(w); };
  auto infWeightIt = std::find_if(allSplitsWeightsSum.begin(), allSplitsWeightsSum.end(), findInf);
  while (infWeightIt != allSplitsWeightsSum.end()) {
    // We have found that the total weight of all the splits for this node is infinite.
    // Therefore, we want to modify the weights to replicate the following
    // sequential behaviors of discrete_distribution_safe/std::discrete_distribution
    // 1) If one or more of the split weights are infinite,
    //    then always pick the first split with infinite weight
    // 2) If none of the split weights are infinite,
    //    then always pick the last split
    auto n = std::distance(allSplitsWeightsSum.begin(), infWeightIt);
    if (myNodeIdx.find(n) != myNodeIdx.end()) {
      LOG_MESSAGE(debug, "Node %u: Handling infinite weights", n);
      static auto nodeSplit = [&n] (const std::pair<uint32_t, double>& split)
                                   { return split.first == n; };
      auto splitIt = std::find_if(mySplitsWeights.begin(), mySplitsWeights.end(), nodeSplit);
      while (splitIt->first == n) {
        // Set all infinite split weights to 1.0 and the finite weights to 0.0
        // This will accomplish the two requirements described above
        splitIt->second = std::isinf(splitIt->second) ? 1.0 : 0.0;
        ++splitIt;
      }
    }
    // Set the total weight for this node to 1.0
    *infWeightIt = 1.0;
    infWeightIt = std::find_if(std::next(infWeightIt), allSplitsWeightsSum.end(), findInf);
  }
  // Then, normalize the weights for all the local splits
  auto normalizeWeight = [&allSplitsWeightsSum] (std::pair<uint32_t, double>& split)
                                                { split.second /= allSplitsWeightsSum[split.first]; };
  std::for_each(mySplitsWeights.begin(), mySplitsWeights.end(), normalizeWeight);
  // Perform a segmented parallel prefix on this vector of <node, normalized weight>
  // with node defining the segment boundaries
  // This will compute the prefix sum of the normalized split weights belonging to every node
  std::vector<std::pair<uint32_t, double>> mySplitsWeightsPrefix(mySplitsWeights.size());
  auto addSplitsWeights = [] (const std::pair<uint32_t, double>& a,
                              const std::pair<uint32_t, double>& b)
                             { return (a.first == b.first) ?
                                      std::make_pair(b.first, a.second + b.second) : b; };
  mxx::global_scan(mySplitsWeights.begin(), mySplitsWeights.end(),
                   mySplitsWeightsPrefix.begin(), addSplitsWeights, false, this->m_comm);
  // Now, we can get the splits for the nodes on this processor
  std::uniform_real_distribution<double> randDist(0.0, 1.0);
  auto first = 0u;
  std::vector<std::tuple<uint32_t, Var, Var, double>> myChosenSplits(myNodeIdx.size() * 2 * numSplits);
  auto chosenIt = myChosenSplits.begin();
  auto g = 0u;
  for (const auto n : myNodeIdx) {
    if (n > g) {
      // Advance the PRNG state to account for the previous nodes
      advance(generator, (n - g) * 2 * numSplits);
      g = n;
    }
    ++g;
    auto nodeSplitsCount = mySplitsCounts[n];
    auto nodeSplitsCountsPrefix = (firstCountsPrefix.first == n) ? firstCountsPrefix.second : 0u;
    auto last = first + nodeSplitsCount - 1;
    auto nodeWeightLower = mySplitsWeightsPrefix[first].second - mySplitsWeights[first].second;
    // There may be some cases in which the last split for a node may
    // not have a cumulative weight of 1.0, because of the following reasons:
    // 1) If all the split weights were finite, but the sum of the weights was infinite
    // 2) Due to eccentricities of floating point arithmetic
    // In such cases, set it to 1.0
    if ((nodeSplitsCountsPrefix + nodeSplitsCount == allSplitsCounts[n]) &&
        std::isless(mySplitsWeightsPrefix[last].second, 1.0)) {
      mySplitsWeightsPrefix[last].second = 1.0;
    }
    auto nodeWeightUpper = mySplitsWeightsPrefix[last].second;
    std::uniform_int_distribution<uint64_t> indexDist(0, allSplitsCounts[n] - 1);
    LOG_MESSAGE(debug, "Node %u: Choosing the splits", n);
    for (auto i = 0u; i < numSplits; ++i) {
      // Pick a split weighted by its score
      auto rand = randDist(generator);
      // Check if the split is local to this processor
      if (std::isgreater(rand, nodeWeightLower) &&
          std::islessequal(rand, nodeWeightUpper)) {
        auto firstSplit = std::next(mySplitsWeightsPrefix.begin(), first);
        auto lastSplit = std::next(firstSplit, nodeSplitsCount);
        auto foundSplit = std::lower_bound(firstSplit, lastSplit, std::make_pair(n, rand));
        auto foundIdx = std::distance(firstSplit, foundSplit);
        auto weightSplit = mySplits[first + foundIdx];
        LOG_MESSAGE(debug, "Chosen parent split using weights: (%s, %g)",
                           this->m_data.varName(std::get<1>(weightSplit)),
                           this->m_data(std::get<1>(weightSplit), std::get<2>(weightSplit)));
        // Index the split with the split index for sorting all the splits for this node later
        *chosenIt = std::make_tuple(i, std::get<1>(weightSplit),
                                    std::get<2>(weightSplit), std::get<3>(weightSplit));
        ++chosenIt;
      }
      else {
        LOG_MESSAGE(debug, "Split for weight %g not in the range [%g, %g)", rand, nodeWeightLower, nodeWeightUpper);
      }
      // Pick a split uniformly at random
      auto randomIdx = indexDist(generator);
      // Check if the split index is local to this processor
      if ((randomIdx >= nodeSplitsCountsPrefix) &&
          (randomIdx < (nodeSplitsCountsPrefix + nodeSplitsCount))) {
        auto randomSplit = mySplits[first + (randomIdx - nodeSplitsCountsPrefix)];
        LOG_MESSAGE(debug, "Chosen parent split at random: (%s, %g)",
                           this->m_data.varName(std::get<1>(randomSplit)),
                           this->m_data(std::get<1>(randomSplit), std::get<2>(randomSplit)));
        // Index the split with the split index for sorting all the splits for this node later
        // Random splits are ordered after all the weight splits
        *chosenIt = std::make_tuple(numSplits + i, std::get<1>(randomSplit),
                                    std::get<2>(randomSplit), std::get<3>(randomSplit));
        ++chosenIt;
      }
      else {
        LOG_MESSAGE(debug, "Random split with index %u not on this processor", randomIdx);
      }
    }
    first += nodeSplitsCount;
  }
  myChosenSplits.resize(std::distance(myChosenSplits.begin(), chosenIt));
  this->m_comm.barrier();
  if (this->m_comm.is_first()) {
    TIMER_ELAPSED("Time taken in choosing splits: ", tChoose);
  }
  TIMER_START(m_tSync);
  // Gather all the chosen splits on all the processors
  // in order to assign them to the corresponding nodes
  // XXX: We can assign the splits for different nodes on different processors
  //      but we are trying to maintain the same state on every processor
  auto allChosenSplits = mxx::allgatherv(myChosenSplits, this->m_comm);
  moduleIt = modules.begin();
  auto countCit = allSplitsCounts.cbegin();
  auto splitIt = allChosenSplits.begin();
  for (auto m = 0u; m < modules.size(); ++m, ++moduleIt) {
    // Synchronize the node parents for this module
    LOG_MESSAGE(info, "Module %u: Synchronizing parents for all nodes", m);
    moduleIt->syncParents(numSplits, countCit, splitIt);
  }
  LOG_MESSAGE(info, "Done synchronizing module parents");
  this->m_comm.barrier();
  TIMER_PAUSE(m_tSync);
}

template <typename Data, typename Var, typename Set>
template <typename Generator>
void
LemonTree<Data, Var, Set>::learnModulesParents_parallel(
  std::list<Module<Data, Var, Set>>& modules,
  Generator& generator,
  const pt::ptree& modulesConfigs
) const
{
  auto regFile = modulesConfigs.get<std::string>("reg_file");
  auto betaMax = modulesConfigs.get<double>("beta_reg");
  auto numSplits = modulesConfigs.get<uint32_t>("num_reg");
  Set candidateParents(this->m_data.numVars());
  if (!regFile.empty()) {
    if (this->m_comm.is_first()) {
      // Read candidate parents from the given file
      this->readCandidateParents(regFile, candidateParents);
    }
    set_bcast(candidateParents, 0, this->m_comm);
  }
  else {
    // Add all the variables as candidate parents
    for (Var v = 0u; v < candidateParents.max(); ++v) {
      candidateParents.insert(v);
    }
  }
  this->learnModulesParents_splits(modules, generator, std::move(candidateParents), betaMax, numSplits);
}

template <typename Data, typename Var, typename Set>
template <typename Generator>
std::list<Module<Data, Var, Set>>
LemonTree<Data, Var, Set>::learnModules(
  const std::multimap<Var, Var>&& coClusters,
  const pt::ptree& modulesConfigs,
  const bool isParallel
) const
{
  auto randomSeed = modulesConfigs.get<uint32_t>("seed");
  Generator generator(randomSeed);
  TIMER_DECLARE(tConstruct);
  auto modules = this->constructModulesWithTrees(std::move(coClusters), generator, modulesConfigs);
  this->m_comm.barrier();
  if (this->m_comm.is_first()) {
    TIMER_ELAPSED("Time taken in learning module trees: ", tConstruct);
  }
  TIMER_DECLARE(tParents);
  if (!isParallel) {
    this->learnModulesParents(modules, generator, modulesConfigs);
  }
  else {
    this->learnModulesParents_parallel(modules, generator, modulesConfigs);
  }
  this->m_comm.barrier();
  if (this->m_comm.is_first()) {
    TIMER_ELAPSED("Time taken in learning module parents: ", tParents);
  }
  return modules;
}

template <typename Data, typename Var, typename Set>
void
LemonTree<Data, Var, Set>::writeParents(
  std::ofstream& stream,
  const std::unordered_map<Var, double>& splits,
  const uint32_t moduleIndex,
  const double cutoff
) const
{
  std::set<std::pair<Var, double>, std::greater<std::pair<Var, double>>> splitsSorted(splits.begin(), splits.end());
  for (const auto& split : splitsSorted) {
    if (split.second > cutoff) {
      stream << this->m_data.varName(split.first) << "\t" << moduleIndex << "\t" << split.second << std::endl;
    }
  }
}

template <typename Data, typename Var, typename Set>
void
LemonTree<Data, Var, Set>::writeModules(
  const std::string& modulesFile,
  const std::list<Module<Data, Var, Set>>& modules
) const
{
  std::string allParentsFile = modulesFile + ".allreg.txt";
  LOG_MESSAGE(info, "Writing all parents to %s", allParentsFile);
  std::ofstream apf(allParentsFile);
  apf.precision(std::numeric_limits<double>::max_digits10);
  std::string topParentsFile = modulesFile + ".topreg.txt";
  LOG_MESSAGE(info, "Writing top 1%% parents to %s", topParentsFile);
  std::ofstream tpf(topParentsFile);
  tpf.precision(std::numeric_limits<double>::max_digits10);
  std::string randParentsFile = modulesFile + ".randomreg.txt";
  LOG_MESSAGE(info, "Writing random parents to %s", randParentsFile);
  std::ofstream rpf(randParentsFile);
  rpf.precision(std::numeric_limits<double>::max_digits10);
  std::string xmlFile = modulesFile + ".xml.gz";
  LOG_MESSAGE(info, "Writing modules to XML file %s", xmlFile);
  std::ofstream file(xmlFile, std::ios_base::out | std::ios_base::binary);
  boost::iostreams::filtering_streambuf<boost::iostreams::output> out;
  out.push(boost::iostreams::gzip_compressor());
  out.push(file);
  std::ostream xmlf(&out);
  xmlf.precision(std::numeric_limits<double>::max_digits10);

  xmlf << "<?xml version='1.0' encoding='iso-8859-1'?>" << std::endl;
  xmlf << "<ModuleNetwork>" << std::endl;
  // First compute the cutoff for top parents
  std::multiset<double, std::greater<double>> allScores;
  for (const auto& module : modules) {
    for (const auto& parent : module.allParents()) {
      allScores.insert(parent.second);
    }
  }
  auto index = static_cast<uint32_t>(round(0.01 * allScores.size()));
  auto cutoff = *std::next(allScores.begin(), index);
  // Now, write the modules
  auto m = 0u;
  for (const auto& module : modules) {
    module.toXML(xmlf, m);
    this->writeParents(apf, module.allParents(), m);
    this->writeParents(tpf, module.allParents(), m, cutoff);
    this->writeParents(rpf, module.randParents(), m);
    ++m;
  }
  xmlf << "</ModuleNetwork>" << std::endl;
}

template <typename Data, typename Var, typename Set>
void
LemonTree<Data, Var, Set>::learnNetwork_sequential(
  const pt::ptree& algoConfigs,
  const std::string& outputDir
) const
{
  using Generator = std::mt19937_64;

  /* GaneSH clustering */
  if (algoConfigs.count("ganesh") == 0) {
    std::cerr << "WARNING: Skipping ganesh and other downstream tasks" << std::endl;
    return;
  }
  TIMER_START(m_tGanesh);
  const auto& ganeshConfigs = algoConfigs.get_child("ganesh");
  auto varClusters = this->clusterVarsGanesh<Generator>(ganeshConfigs);
  TIMER_PAUSE(m_tGanesh);
  auto clusterFile = ganeshConfigs.get<std::string>("output_file");
  if (!clusterFile.empty()) {
    TIMER_START(m_tWrite);
    clusterFile = outputDir + "/" + clusterFile;
    this->writeVarClusters(clusterFile, varClusters);
    TIMER_PAUSE(m_tWrite);
  }

  /* Consensus clustering */
  if (algoConfigs.count("tight_clusters") == 0) {
    std::cerr << "WARNING: Skipping tight_clusters and other downstream tasks" << std::endl;
    return;
  }
  TIMER_START(m_tConsensus);
  const auto& consensusConfigs = algoConfigs.get_child("tight_clusters");
  auto coClusters = this->clusterConsensus(std::move(varClusters), consensusConfigs);
  TIMER_PAUSE(m_tConsensus);
  auto consensusFile = consensusConfigs.get<std::string>("output_file");
  if (!consensusFile.empty()) {
    TIMER_START(m_tWrite);
    consensusFile = outputDir + "/" + consensusFile;
    this->writeConsensusCluster(consensusFile, coClusters);
    TIMER_PAUSE(m_tWrite);
  }

  /* Module learning */
  if (algoConfigs.count("regulators") == 0) {
    std::cerr << "WARNING: Skipping regulators and other downstream tasks" << std::endl;
    return;
  }
  TIMER_START(m_tModules);
  const auto& modulesConfigs = algoConfigs.get_child("regulators");
  auto modules = this->learnModules<Generator>(std::move(coClusters), modulesConfigs);
  TIMER_PAUSE(m_tModules);
  auto modulesFile = modulesConfigs.get<std::string>("output_file");
  if (!modulesFile.empty()) {
    TIMER_START(m_tWrite);
    modulesFile = outputDir + "/" + modulesFile;
    this->writeModules(modulesFile, modules);
    TIMER_PAUSE(m_tWrite);
  }
}

template <typename Data, typename Var, typename Set>
void
LemonTree<Data, Var, Set>::learnNetwork_parallel(
  const pt::ptree& algoConfigs,
  const std::string& outputDir
) const
{
  using Generator = std::mt19937_64;

  /* GaneSH clustering */
  if (algoConfigs.count("ganesh") == 0) {
    if (this->m_comm.is_first()) {
      std::cerr << "WARNING: Skipping ganesh and other downstream tasks" << std::endl;
    }
    return;
  }
  TIMER_START(m_tGanesh);
  const auto& ganeshConfigs = algoConfigs.get_child("ganesh");
  auto varClusters = this->clusterVarsGanesh<Generator>(ganeshConfigs);
  this->m_comm.barrier();
  TIMER_PAUSE(m_tGanesh);
  auto clusterFile = ganeshConfigs.get<std::string>("output_file");
  if (!clusterFile.empty() && this->m_comm.is_first()) {
    TIMER_START(m_tWrite);
    clusterFile = outputDir + "/" + clusterFile;
    this->writeVarClusters(clusterFile, varClusters);
    TIMER_PAUSE(m_tWrite);
  }
  this->m_comm.barrier();

  /* Consensus clustering */
  if (algoConfigs.count("tight_clusters") == 0) {
    if (this->m_comm.is_first()) {
      std::cerr << "WARNING: Skipping tight_clusters and other downstream tasks" << std::endl;
    }
    return;
  }
  TIMER_START(m_tConsensus);
  const auto& consensusConfigs = algoConfigs.get_child("tight_clusters");
  auto coClusters = this->clusterConsensus(std::move(varClusters), consensusConfigs);
  TIMER_PAUSE(m_tConsensus);
  auto consensusFile = consensusConfigs.get<std::string>("output_file");
  if (!consensusFile.empty() && this->m_comm.is_first()) {
    TIMER_START(m_tWrite);
    consensusFile = outputDir + "/" + consensusFile;
    this->writeConsensusCluster(consensusFile, coClusters);
    TIMER_PAUSE(m_tWrite);
  }
  this->m_comm.barrier();

  /* Module learning */
  if (algoConfigs.count("regulators") == 0) {
    if (this->m_comm.is_first()) {
      std::cerr << "WARNING: Skipping regulators and other downstream tasks" << std::endl;
    }
    return;
  }
  TIMER_START(m_tModules);
  const auto& modulesConfigs = algoConfigs.get_child("regulators");
  auto modules = this->learnModules<Generator>(std::move(coClusters), modulesConfigs, true);
  this->m_comm.barrier();
  TIMER_PAUSE(m_tModules);
  auto modulesFile = modulesConfigs.get<std::string>("output_file");
  if (!modulesFile.empty() && this->m_comm.is_first()) {
    TIMER_START(m_tWrite);
    modulesFile = outputDir + "/" + modulesFile;
    this->writeModules(modulesFile, modules);
    TIMER_PAUSE(m_tWrite);
  }
}

#endif // DETAIL_LEMONTREE_HPP_
