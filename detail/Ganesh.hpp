/**
 * @file Ganesh.hpp
 * @brief Implementation of functionality for GaneSH clustering.
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
#ifndef DETAIL_GANESH_HPP_
#define DETAIL_GANESH_HPP_

#include "PrimaryCluster.hpp"


/**
 * @brief Class that implements the two-way Gibbs clustering algorithm,
 *        called GaneSH, by Joshi et al.
 *
 * @tparam Data Type of the object which is used for querying the data.
 * @tparam Var Type of variable indices (expected to be an integer type).
 * @tparam Set Type of set container.
 */
template <typename Data, typename Var, typename Set>
class Ganesh {
public:
  Ganesh(const Data&);

  ~Ganesh();

  template <typename Generator>
  void
  initializeRandom(Generator&, const Var);

  template <typename Generator>
  void
  initializeGiven(Generator&, const std::list<Set>&);

  template <typename Generator>
  void
  clusterTwoWay(Generator&, const mxx::comm&);

  template <typename Generator>
  void
  clusterSecondary(Generator&, const mxx::comm* const = nullptr, const uint32_t = 1);

  const std::list<PrimaryCluster<Data, Var, Set>>&
  primaryClusters() const;

private:
  template <typename Generator>
  void
  initializeSecondaryRandom(Generator&, const Var);

  void
  removeEmptyClusters();

  template <typename Generator>
  Var
  chooseReassignCluster(Generator&, const Var, const double);

  template <typename Generator>
  Var
  chooseReassignCluster(Generator&, const mxx::comm&, const Var, const double);

  template <typename Generator>
  void
  reassignPrimary(Generator&, const mxx::comm&, const Var);

  template <typename Generator>
  Var
  chooseMergeCluster(Generator&, const typename std::list<PrimaryCluster<Data, Var, Set>>::iterator&);

  template <typename Generator>
  Var
  chooseMergeCluster(Generator&, const mxx::comm&, const typename std::list<PrimaryCluster<Data, Var, Set>>::iterator&);

  template <typename Generator>
  bool
  mergeCluster(Generator&, const mxx::comm&, typename std::list<PrimaryCluster<Data, Var, Set>>::iterator&);

  template <typename Generator>
  void
  clusterPrimary(Generator&, const mxx::comm&);

private:
  std::list<PrimaryCluster<Data, Var, Set>> m_cluster;
  std::vector<typename std::list<PrimaryCluster<Data, Var, Set>>::iterator> m_membership;
  const Data& m_data;
}; // class Ganesh

template <typename Data, typename Var, typename Set>
/**
 * @brief Constructs a Gibbs clustering object.
 *
 * @param data The data provider.
 */
Ganesh<Data, Var, Set>::Ganesh(
  const Data& data
) : m_cluster(),
    m_membership(),
    m_data(data)
{
}

template <typename Data, typename Var, typename Set>
/**
 * @brief Default destructor.
 */
Ganesh<Data, Var, Set>::~Ganesh(
)
{
}

template <typename Data, typename Var, typename Set>
/**
 * @brief Randomly initializes the secondary clusters
 *        for all the primary clusters.
 *
 * @tparam Generator Type of PRNG used for generating random numbers.
 * @param generator Reference to the instance of the PRNG.
 */
template <typename Generator>
void
Ganesh<Data, Var, Set>::initializeSecondaryRandom(
  Generator& generator,
  const Var numSecondaryClusters
)
{
  for (auto& cluster : m_cluster) {
    cluster.randomSecondary(generator, numSecondaryClusters);
  }
}

template <typename Data, typename Var, typename Set>
/**
 * @brief Randomly initializes the primary clusters and
 *        the corresponding secondary clusters.
 *
 * @tparam Generator Type of PRNG used for generating random numbers.
 * @param generator Reference to the instance of the PRNG.
 */
template <typename Generator>
void
Ganesh<Data, Var, Set>::initializeRandom(
  Generator& generator,
  const Var numPrimary
)
{
  const auto n = m_data.numVars();
  const auto m = m_data.numObs();
  LOG_MESSAGE(info, "Randomly assigning primary variables to %u clusters", static_cast<uint32_t>(numPrimary));
  std::vector<PrimaryCluster<Data, Var, Set>> cluster(numPrimary, PrimaryCluster<Data, Var, Set>(m_data, n, m));
  std::uniform_int_distribution<Var> clusterDistrib(0, numPrimary - 1);
  for (Var e = 0; e < n; ++e) {
    // Pick a cluster uniformly at random
    auto c = clusterDistrib(generator);
    cluster[c].insert(e);
  }
  m_cluster = std::list<PrimaryCluster<Data, Var, Set>>(cluster.begin(), cluster.end());
  this->removeEmptyClusters();
  m_membership.resize(n, m_cluster.end());
  for (auto cIt = m_cluster.begin(); cIt != m_cluster.end(); ++cIt) {
    for (const auto e : cIt->elements()) {
      m_membership[e] = cIt;
    }
  }
  LOG_MESSAGE(info, "Assigned primary variables to %u clusters", m_cluster.size());
  auto numSecondaryClusters = static_cast<Var>(sqrt(m));
  this->initializeSecondaryRandom(generator, numSecondaryClusters);
}

template <typename Data, typename Var, typename Set>
/**
 * @brief Initializes the primary clusters with the given clusters
 *        and randomly initializes the corresponding secondary clusters.
 */
template <typename Generator>
void
Ganesh<Data, Var, Set>::initializeGiven(
  Generator& generator,
  const std::list<Set>& givenClusters
)
{
  const auto m = m_data.numObs();
  LOG_MESSAGE(info, "Randomly initializing using the given primary clusters");
  m_cluster.clear();
  for (const auto& cluster : givenClusters) {
    m_cluster.emplace_back(m_data, cluster, m);
  }
  auto numSecondaryClusters = static_cast<Var>(sqrt(m));
  this->initializeSecondaryRandom(generator, numSecondaryClusters);
}

template <typename Data, typename Var, typename Set>
/**
 * @brief Removes all the empty primary clusters.
 */
void
Ganesh<Data, Var, Set>::removeEmptyClusters(
)
{
  auto emptyCluster = [] (const PrimaryCluster<Data, Var, Set>& cluster)
                         { return cluster.empty(); };
  m_cluster.remove_if(emptyCluster);
}

template <typename Data, typename Var, typename Set>
template <typename Generator>
Var
Ganesh<Data, Var, Set>::chooseReassignCluster(
  Generator& generator,
  const Var given,
  const double singleScore
)
{
  // Compute the weight of the var being in its separate cluster
  // as well as all the existing clusters
  std::vector<double> weight(m_cluster.size() + 1);
  weight[0] = 1.0;
  auto maxDiff = weight[0];
  auto wIt = weight.begin() + 1;
  // Only compute score diffs for existing clusters
  for (auto cIt = m_cluster.begin(); cIt != m_cluster.end(); ++cIt, ++wIt) {
    auto thisDiff = cIt->scoreInsertPrimary(given) -
                    (cIt->score() + singleScore);
    *wIt = thisDiff;
    maxDiff = std::max(thisDiff, maxDiff);
  }
  for (auto& w : weight) {
    w = exp(w - maxDiff);
  }
  // Pick a cluster using the computed weights
  auto distrib = discrete_distribution_safe<Var>(weight.cbegin(), weight.cend());
  return distrib(generator);
}

template <typename Data, typename Var, typename Set>
template <typename Generator>
Var
Ganesh<Data, Var, Set>::chooseReassignCluster(
  Generator& generator,
  const mxx::comm& comm,
  const Var given,
  const double singleScore
)
{
  mxx::blk_dist block(m_cluster.size(), comm.size(), comm.rank());
  // Compute the weight of the var being in its separate cluster
  // as well as all the existing clusters
  std::vector<double> myWeights(block.local_size() + static_cast<uint8_t>(comm.is_first()));
  auto myMaxDiff = std::numeric_limits<double>::lowest();
  if (comm.is_first()) {
    myWeights[0] = 1.0;
    myMaxDiff = 1.0;
  }
  auto wIt = std::next(myWeights.begin(), static_cast<uint8_t>(comm.is_first()));
  // Only compute score diffs for existing clusters
  auto cIt = std::next(m_cluster.begin(), block.eprefix_size());
  for (auto c = block.eprefix_size(); c < block.iprefix_size(); ++c, ++cIt, ++wIt) {
    auto thisDiff = cIt->scoreInsertPrimary(given) -
                    (cIt->score() + singleScore);
    *wIt = thisDiff;
    myMaxDiff = std::max(thisDiff, myMaxDiff);
  }
  auto allMaxDiff = mxx::allreduce(myMaxDiff, mxx::max<double>(), comm);
  for (auto& w : myWeights) {
    w = exp(w - allMaxDiff);
  }
  auto allWeights = mxx::allgatherv(myWeights, comm);
  // Pick a cluster using the computed weights
  auto distrib = discrete_distribution_safe<Var>(allWeights.cbegin(), allWeights.cend());
  return distrib(generator);
}

template <typename Data, typename Var, typename Set>
/**
 * @brief Moves the given primary variable to a different
 *        primary cluster.
 *
 * @tparam Generator Type of PRNG used for generating random numbers.
 * @param generator Reference to the instance of the PRNG.
 * @param given The index of the primary variable to be moved.
 */
template <typename Generator>
void
Ganesh<Data, Var, Set>::reassignPrimary(
  Generator& generator,
  const mxx::comm& comm,
  const Var given
)
{
  LOG_MESSAGE(debug, "Reassigning primary variable %u", static_cast<uint32_t>(given));
  auto oldCluster = m_membership[given];
  m_membership[given] = m_cluster.end();
  // Create a copy of the old cluster to get the same
  // clustering of the secondary elements
  PrimaryCluster<Data, Var, Set> newCluster(*oldCluster);
  // Remove all the primary elements from the copy
  newCluster.scoreClear();
  newCluster.clear();
  // Insert this element as the only primary member of the copy
  newCluster.insert(given);
  if (oldCluster->size() > 1) {
    // Remove the element and update the score of the cluster
    oldCluster->scoreErasePrimary(given, true);
    oldCluster->erase(given);
  }
  else {
    // Remove the cluster if the var was its only element
    LOG_MESSAGE(debug, "Removing the old cluster of the variable");
    m_cluster.erase(oldCluster);
  }
  auto c = m_cluster.size() + 1;
  if (comm.size() == 1) {
    c = this->chooseReassignCluster(generator, given, newCluster.score());
  }
  else {
    c = this->chooseReassignCluster(generator, comm, given, newCluster.score());
  }
  if (c == 0) {
    // The variable will stay in its own cluster
    LOG_MESSAGE(info, "Primary variable %u assigned to a newly created cluster", static_cast<uint32_t>(given));
    newCluster.clearSecondary();
    newCluster.randomSecondary(generator, static_cast<Var>(sqrt(m_data.numObs())));
    m_cluster.push_back(std::move(newCluster));
    m_membership[given] = std::prev(m_cluster.end());
  }
  else {
    // Add the variable to the chosen cluster
    LOG_MESSAGE(info, "Primary variable %u assigned to the existing cluster %u",
                      static_cast<uint32_t>(given), static_cast<uint32_t>(c - 1));
    auto chosen = std::next(m_cluster.begin(), c - 1);
    chosen->scoreInsertPrimary(given, true);
    chosen->insert(given);
    m_membership[given] = chosen;
  }
}

template <typename Data, typename Var, typename Set>
template <typename Generator>
Var
Ganesh<Data, Var, Set>::chooseMergeCluster(
  Generator& generator,
  const typename std::list<PrimaryCluster<Data, Var, Set>>::iterator& given
)
{
  auto givenScore = given->score();
  // Compute the weight of merging this cluster with
  // all the other clusters which are not empty
  std::vector<double> weight(m_cluster.size(), 0.0);
  auto wIt = weight.begin();
  for (auto cIt = m_cluster.begin(); cIt != m_cluster.end(); ++cIt, ++wIt) {
    if (cIt != given) {
      PrimaryCluster<Data, Var, Set> merged(*given, *cIt);
      auto thisDiff = merged.score() - (cIt->score() + givenScore);
      *wIt = exp(thisDiff);
    }
    else {
      *wIt = 1.0;
    }
  }
  // Choose a cluster using the computed weights
  auto distrib = discrete_distribution_safe<Var>(weight.cbegin(), weight.cend());
  return distrib(generator);
}

template <typename Data, typename Var, typename Set>
template <typename Generator>
Var
Ganesh<Data, Var, Set>::chooseMergeCluster(
  Generator& generator,
  const mxx::comm& comm,
  const typename std::list<PrimaryCluster<Data, Var, Set>>::iterator& given
)
{
  auto givenScore = given->score();
  mxx::blk_dist block(m_cluster.size(), comm.size(), comm.rank());
  std::vector<double> myWeights(block.local_size());
  auto wIt = myWeights.begin();
  auto cIt = std::next(m_cluster.begin(), block.eprefix_size());
  for (auto c = block.eprefix_size(); c < block.iprefix_size(); ++c, ++cIt, ++wIt) {
    if (cIt != given) {
      PrimaryCluster<Data, Var, Set> merged(*given, *cIt);
      auto thisDiff = merged.score() - (cIt->score() + givenScore);
      *wIt = exp(thisDiff);
    }
    else {
      *wIt = 1.0;
    }
  }
  auto allWeights = mxx::allgatherv(myWeights, comm);
  // Choose a cluster using all the computed weights
  auto distrib = discrete_distribution_safe<Var>(allWeights.cbegin(), allWeights.cend());
  return distrib(generator);
}

template <typename Data, typename Var, typename Set>
/**
 * @brief Merges the given primary cluster with another primary cluster.
 *
 * @tparam Generator Type of PRNG used for generating random numbers.
 * @param generator Reference to the instance of the PRNG.
 */
template <typename Generator>
bool
Ganesh<Data, Var, Set>::mergeCluster(
  Generator& generator,
  const mxx::comm& comm,
  typename std::list<PrimaryCluster<Data, Var, Set>>::iterator& given
)
{
  auto c = m_cluster.size();
  if (comm.size() == 1) {
    c = this->chooseMergeCluster(generator, given);
  }
  else {
    c = this->chooseMergeCluster(generator, comm, given);
  }
  auto chosen = std::next(m_cluster.begin(), c);
  if (chosen != given) {
    LOG_MESSAGE(info, "Merging given cluster with cluster %u", static_cast<uint32_t>(c));
    // Merge this cluster with the chosen cluster
    // and update the membership of all the moved elements
    for (const auto e : given->elements()) {
      m_membership[e] = chosen;
    }
    chosen->scoreMerge(*given, true);
    chosen->merge(*given);
    return true;
  }
  else {
    LOG_MESSAGE(info, "Not merging given cluster");
    return false;
  }
}

template <typename Data, typename Var, typename Set>
/**
 * @brief Performs a Gibbs clustering step for the primary variables.
 *
 * @tparam Generator Type of PRNG used for generating random numbers.
 * @param generator Reference to the instance of the PRNG.
 */
template <typename Generator>
void
Ganesh<Data, Var, Set>::clusterPrimary(
  Generator& generator,
  const mxx::comm& comm
)
{
  const auto n = m_data.numVars();
  std::uniform_int_distribution<Var> varDistrib(0, n - 1);
  // Reassign a random variable for n iterations
  LOG_MESSAGE(info, "Reassigning primary variables");
  for (auto i = 0u; i < n; ++i) {
    auto v = varDistrib(generator);
    this->reassignPrimary(generator, comm, v);
  }
  LOG_MESSAGE(info, "Done reassigning primary variables");
  // Try to merge clusters
  LOG_MESSAGE(info, "Merging primary clusters (number of clusters = %u)", m_cluster.size());
  for (auto cIt = m_cluster.begin(); (cIt != m_cluster.end()) && (m_cluster.size() > 1); ) {
    if (this->mergeCluster(generator, comm, cIt)) {
      cIt = m_cluster.erase(cIt);
    }
    else {
      ++cIt;
    }
  }
  LOG_MESSAGE(info, "Done merging primary clusters (number of clusters = %u)", m_cluster.size());
}

template <typename Data, typename Var, typename Set>
/**
 * @brief Performs a Gibbs clustering step for both primary as
 *        well as secondary variables.
 *
 * @tparam Generator Type of PRNG used for generating random numbers.
 * @param generator Reference to the instance of the PRNG.
 */
template <typename Generator>
void
Ganesh<Data, Var, Set>::clusterTwoWay(
  Generator& generator,
  const mxx::comm& comm
)
{
  LOG_MESSAGE(info, "Clustering primary variables");
  this->clusterPrimary(generator, comm);
  LOG_MESSAGE(info, "Done clustering primary variables");
  this->clusterSecondary(generator, &comm, 50);
}

template <typename Data, typename Var, typename Set>
/**
 * @brief Performs a Gibbs clustering step for secondary variables.
 *
 * @tparam Generator Type of PRNG used for generating random numbers.
 * @param generator Reference to the instance of the PRNG.
 * @param numReps Number of times clustering of secondary variables
 *                should be repeated.
 */
template <typename Generator>
void
Ganesh<Data, Var, Set>::clusterSecondary(
  Generator& generator,
  const mxx::comm* const comm,
  const uint32_t numReps
)
{
  LOG_MESSAGE(info, "Clustering secondary variables for all the primary clusters");
  if ((comm != nullptr) && (comm->size() > 1) && (m_cluster.size() > 1)) {
    auto perClusterGenerated = m_data.numObs() * numReps * 3;
    if (static_cast<uint32_t>(comm->size()) > m_cluster.size()) {
      LOG_MESSAGE(debug, "Clustering in parallel by splitting communicator");
      // Split the communicator with more than one process per cluster
      mxx::blk_dist commBlock(comm->size(), m_cluster.size(), 0);
      auto myCluster = commBlock.rank_of(comm->rank());
      auto clusterComm = comm->split(myCluster);
      // Advance the generator state for previous clusters
      ::advance(generator, myCluster * perClusterGenerated);
      auto cIt = std::next(m_cluster.begin(), myCluster);
      // First, learn secondary clusters for the local primary cluster
      // Each process will call clusterSecondary for just one cluster
      cIt->clusterSecondary(generator, &clusterComm, numReps);
      // Then, synchronize secondary clusters for all the primary clusters
      cIt = m_cluster.begin();
      for (auto c = 0u; c < m_cluster.size(); ++c, ++cIt) {
        // Synchronize the cluster using the first process for the
        // cluster as the source
        cIt->syncSecondary(*comm, commBlock.eprefix_size(c));
      }
      // Advance the generator state for next clusters
      ::advance(generator, (m_cluster.size() - myCluster - 1) * perClusterGenerated);
    }
    else {
      LOG_MESSAGE(debug, "Clustering in parallel by splitting clusters");
      // Assign one or more clusters per process; no need to split the communicator
      mxx::blk_dist clusterBlock(m_cluster.size(), comm->size(), comm->rank());
      // Advance the generator state for previous clusters
      ::advance(generator, clusterBlock.eprefix_size() * perClusterGenerated);
      // First, learn secondary clusters for all the local primary clusters
      auto cIt = std::next(m_cluster.begin(), clusterBlock.eprefix_size());
      for (auto c = clusterBlock.eprefix_size(); c < clusterBlock.iprefix_size(); ++c, ++cIt) {
        cIt->clusterSecondary(generator, nullptr, numReps);
      }
      // Advance the generator state for next clusters
      ::advance(generator, (m_cluster.size() - clusterBlock.iprefix_size()) * perClusterGenerated);
      // Then, synchronize secondary clusters for all the primary clusters
      cIt = m_cluster.begin();
      for (auto c = 0u; c < m_cluster.size(); ++c, ++cIt) {
        cIt->syncSecondary(*comm, clusterBlock.rank_of(c));
      }
    }
  }
  else {
    // There is no way to learn different secondary clusters in parallel
    // We may compute the scores for reassignments and merges in parallel
    for (auto& cluster : m_cluster) {
      cluster.clusterSecondary(generator, comm, numReps);
    }
  }
  LOG_MESSAGE(info, "Done clustering secondary variables");
}

template <typename Data, typename Var, typename Set>
/**
 * @brief Returns the primary clusters.
 */
const std::list<PrimaryCluster<Data, Var, Set>>&
Ganesh<Data, Var, Set>::primaryClusters(
) const
{
  return m_cluster;
}

#endif // DETAIL_GANESH_HPP_
