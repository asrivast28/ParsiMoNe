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

#include "utils/Random.hpp"


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
  clusterTwoWay(Generator&);

  template <typename Generator>
  void
  clusterSecondary(Generator&, const uint32_t = 1);

  const std::list<PrimaryCluster<Data, Var, Set>>&
  primaryClusters() const;

private:
  template <typename Generator>
  void
  initializeSecondaryRandom(Generator&, const Var);

  void
  removeEmptyClusters();

  template <typename Generator>
  void
  reassignPrimary(Generator&, const Var);

  double
  scoreMerge(PrimaryCluster<Data, Var, Set>&, PrimaryCluster<Data, Var, Set>&);

  template <typename Generator>
  bool
  mergeCluster(Generator&, typename std::list<PrimaryCluster<Data, Var, Set>>::iterator&);

  template <typename Generator>
  void
  clusterPrimary(Generator&);

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
  for (auto cit = m_cluster.begin(); cit != m_cluster.end(); ++cit) {
    for (const auto e : cit->elements()) {
      m_membership[e] = cit;
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
  // Compute the weight of the var being in its separate cluster
  // as well as all the existing clusters
  std::vector<double> weight(m_cluster.size() + 1);
  weight[0] = 1.0;
  auto maxDiff = weight[0];
  auto singleScore = newCluster.score();
  auto wit = weight.begin() + 1;
  // Only compute score diffs for existing clusters
  for (auto& cluster : m_cluster) {
    auto thisDiff = cluster.scoreInsertPrimary(given) -
                    cluster.score() - singleScore;
    *wit = thisDiff;
    ++wit;
    maxDiff = std::max(thisDiff, maxDiff);
  }
  for (auto& w : weight) {
    w = exp(w - maxDiff);
  }
  // Pick a cluster using the computed weights
  auto c = discrete_distribution_pick<Var>(weight.cbegin(), weight.cend(), generator);
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
/**
 * @brief Computes the score of merging two primary clusters.
 *        The merged primary cluster has just one secondary cluster.
 *
 * @param first The first primary cluster to be merged.
 * @param second The second primary cluster to be merged.
 *
 * @return The score for merging the two given primary clusters.
 */
double
Ganesh<Data, Var, Set>::scoreMerge(
  PrimaryCluster<Data, Var, Set>& first,
  PrimaryCluster<Data, Var, Set>& second
)
{
  auto firstScore = first.score();
  auto secondScore = second.score();
  PrimaryCluster<Data, Var, Set> merged(first, second);
  auto mergeScore = merged.score();
  auto scoreDiff = mergeScore - (firstScore + secondScore);
  return scoreDiff;
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
  typename std::list<PrimaryCluster<Data, Var, Set>>::iterator& given
)
{
  // Compute the weight of merging this cluster with
  // all the other clusters which are not empty
  std::vector<double> weight(m_cluster.size(), 0.0);
  auto wit = weight.begin();
  for (auto cit = m_cluster.begin(); cit != m_cluster.end(); ++cit, ++wit) {
    if (cit != given) {
      *wit = exp(this->scoreMerge(*given, *cit));
    }
    else {
      *wit = 1.0;
    }
  }

  // Choose a cluster using the computed weights
  auto c = discrete_distribution_pick<Var>(weight.cbegin(), weight.cend(), generator);
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
  Generator& generator
)
{
  const auto n = m_data.numVars();
  std::uniform_int_distribution<Var> varDistrib(0, n - 1);
  // Reassign a random variable for n iterations
  LOG_MESSAGE(info, "Reassigning primary variables");
  for (auto i = 0u; i < n; ++i) {
    auto v = varDistrib(generator);
    this->reassignPrimary(generator, v);
  }
  LOG_MESSAGE(info, "Done reassigning primary variables");
  // Try to merge clusters
  LOG_MESSAGE(info, "Merging primary clusters (number of clusters = %u)", m_cluster.size());
  for (auto cit = m_cluster.begin(); (cit != m_cluster.end()) && (m_cluster.size() > 1); ) {
    if (this->mergeCluster(generator, cit)) {
      cit = m_cluster.erase(cit);
    }
    else {
      ++cit;
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
  Generator& generator
)
{
  LOG_MESSAGE(info, "Clustering primary variables");
  this->clusterPrimary(generator);
  LOG_MESSAGE(info, "Done clustering primary variables");
  this->clusterSecondary(generator, 50);
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
  const uint32_t numReps
)
{
  LOG_MESSAGE(info, "Clustering secondary variables for all the primary clusters");
  for (auto& cluster : m_cluster) {
    cluster.clusterSecondary(generator, numReps);
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
