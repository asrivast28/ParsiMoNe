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
#include "Generator.hpp"

#include "mxx/partition.hpp"

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
std::list<std::list<Set>>
LemonTree<Data, Var, Set>::clusterVarsGanesh(
  const pt::ptree& ganeshConfigs
) const
{
  auto randomSeed = ganeshConfigs.get<uint32_t>("seed");
  auto initClusters = ganeshConfigs.get<Var>("init_num_clust");
  if ((initClusters == 0) || (initClusters > this->m_data.numVars())) {
    initClusters = this->m_data.numVars() / 2;
  }
  auto numRuns = ganeshConfigs.get<uint32_t>("num_runs");
  auto numSteps = ganeshConfigs.get<uint32_t>("num_steps");
  std::list<std::list<Set>> sampledClusters;
  Ganesh<Data, Var, Set> ganesh(this->m_data);
  for (auto r = 0u; r < numRuns; ++r) {
    LOG_MESSAGE(info, "Run %u", r);
    // XXX: Seeding PRNG in this way to compare the results
    //      with Lemon-Tree; the seeding can be moved outside later
    Generator generator(randomSeed + r);
    ganesh.initializeRandom(generator, initClusters);
    for (auto s = 0u; s <= numSteps; ++s) {
      LOG_MESSAGE(info, "Step %u", s);
      ganesh.clusterTwoWay(generator);
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
  OptimalBeta ob(0.0, betaMax, 1e-5);
  std::vector<uint32_t> moduleNodeCount(modules.size());
  auto totalNodes = 0u;
  auto moduleCit = modules.cbegin();
  for (auto m = 0u; m < modules.size(); ++m, ++moduleCit) {
    moduleNodeCount[m] = moduleCit->nodeCount();
    totalNodes += moduleNodeCount[m];
  }
  std::vector<uint32_t> moduleNodeCountPrefix(moduleNodeCount.size());
  std::partial_sum(moduleNodeCount.cbegin(), moduleNodeCount.cend(), moduleNodeCountPrefix.begin());
  mxx::blk_dist block(totalNodes, this->m_comm.size(), this->m_comm.rank());
  // First, advance the PRNG state to account for the number of
  // random generations for nodes on the previous ranks
  advance(generator, block.eprefix_size() * 2 * numSplits);
  // XXX: We need to track the validity of nodes because some nodes may not learn any splits
  //      and that information will be required for synchornization later
  std::vector<uint8_t> myValidNodes(block.local_size());
  std::vector<std::tuple<Var, Var, double>> myNodeSplits(block.local_size() * 2 * numSplits);
  // Find the indices for the modules which contain the first and last node on this rank
  auto myFirst = std::distance(moduleNodeCountPrefix.cbegin(),
                               std::lower_bound(moduleNodeCountPrefix.cbegin(),
                                                moduleNodeCountPrefix.cend(),
                                                block.eprefix_size() + 1));
  auto myLast = std::distance(moduleNodeCountPrefix.cbegin(),
                              std::lower_bound(moduleNodeCountPrefix.cbegin(),
                                               moduleNodeCountPrefix.cend(),
                                               block.iprefix_size()));
  auto moduleIt = std::next(modules.begin(), myFirst);
  auto validIt = myValidNodes.begin();
  auto splitIt = myNodeSplits.begin();
  auto prevNodeCount = block.eprefix_size();
  for (auto m = myFirst; m <= myLast; ++m, ++moduleIt) {
    auto firstNode = prevNodeCount - (moduleNodeCountPrefix[m] - moduleNodeCount[m]);
    auto numNodes = std::min(static_cast<uint32_t>(block.iprefix_size()), moduleNodeCountPrefix[m]) - prevNodeCount;
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
