/**
 * @file TreeNode.hpp
 * @brief Implementation of the module tree node class.
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
#ifndef DETAIL_TREENODE_HPP_
#define DETAIL_TREENODE_HPP_

#include "Assignment.hpp"
#include "Random.hpp"

#include <trng/uniform_int_dist.hpp>


class OptimalBeta {
public:
  OptimalBeta(
    const double min,
    const double max,
    const double acc
  ) : m_min(min),
      m_max(max),
      m_acc(acc)
  {
  }

  template <typename Data, typename Var, typename Set>
  double
  find(
    const Assignment<Data, Var, Set>& assmt,
    const double sv,
    const int sign
  ) const
  {
    auto min = m_min;
    auto max = m_max;
    auto optimal = 0.0;
    auto found = this->init(assmt, sv, sign, min, max);
    if (found) {
      found = this->bisect(assmt, sv, sign, min, max, optimal);
    }
    return found ? optimal : std::nan("");
  }

private:
  template <typename Data, typename Var, typename Set>
  bool
  init(
    const Assignment<Data, Var, Set>& assmt,
    const double sv,
    const int sign,
    double& min,
    double& max
  ) const
  {
    auto f1 = assmt.evaluate(sv, sign, min);
    auto f2 = assmt.evaluate(sv, sign, max);
    auto found = false;
    for (auto i = 0u; i < m_initTries; ++i) {
      if (std::isless(f1 * f2, 0)) {
        found = true;
        break;
      }
      min = max;
      max = m_factor * max;
      f2 = assmt.evaluate(sv, sign, max);
    }
    LOG_MESSAGE_IF(!found, trace, "Unable to initialize beta for split (%s, %g)",
                                  assmt.varName(), sv);
    return found;
  }

  template <typename Data, typename Var, typename Set>
  bool
  bisect(
    const Assignment<Data, Var, Set>& assmt,
    const double sv,
    const int sign,
    const double min,
    const double max,
    double& optimal
  ) const
  {
    auto fMin = assmt.evaluate(sv, sign, min);
    auto diff = 0.0;
    if (std::isless(fMin, 0)) {
      optimal = min;
      diff = max - min;
    }
    else {
      optimal = max;
      diff = min - max;
    }
    auto found = false;
    for (auto i = 0u; i < m_bisectTries; ++i) {
      diff *= 0.5;
      auto mid = optimal + diff;
      auto fMid = assmt.evaluate(sv, sign, mid);
      if (std::islessequal(fMid, 0)) {
        optimal = mid;
      }
      if (std::isless(abs(diff), m_acc) || (fMid == 0)) {
        found = true;
        break;
      }
    }
    LOG_MESSAGE_IF(!found, warning, "Unable to find beta for split (%s, %g)",
                                    assmt.varName(), sv);
    return found;
  }

private:
  static constexpr double m_factor = 2.0;
  static constexpr uint32_t m_initTries = 50;
  static constexpr uint32_t m_bisectTries = 40;

private:
  const double m_min;
  const double m_max;
  const double m_acc;
}; // class OptimalBeta


template <typename Data, typename Var, typename Set>
class TreeNode {
public:
  TreeNode(const Data&, const Set&, const Set&);

  TreeNode(const std::shared_ptr<TreeNode<Data, Var, Set>>&, const std::shared_ptr<TreeNode<Data, Var, Set>>&);

  const std::pair<std::shared_ptr<TreeNode<Data, Var, Set>>, std::shared_ptr<TreeNode<Data, Var, Set>>>&
  children() const;

  void
  makeLeaf();

  uint32_t
  nodeCount(const bool = false) const;

  std::list<TreeNode<Data, Var, Set>*>
  nodes(const bool = false);

  const Set&
  observations() const;

  double
  score() const;

  double
  mergeScore(const bool) const;

  double
  mean() const;

  void
  prune(const double);

  uint64_t
  splitWeight(const Set&) const;

  std::list<std::tuple<Var, Var, double>>
  candidateParentsSplits(const Set&, const OptimalBeta&, const uint64_t, const uint64_t) const;

  template <typename Generator>
  bool
  learnParentsSplits(Generator&, const Set&, const OptimalBeta&, const uint32_t);

  template <typename Generator, typename SplitIt>
  bool
  learnParentsSplits(Generator&, const Set&, const OptimalBeta&, const uint32_t, SplitIt, SplitIt) const;

  template <typename SplitIt>
  void
  setWeightSplits(const SplitIt&, const SplitIt&);

  template <typename SplitIt>
  void
  setRandomSplits(const SplitIt&, const SplitIt&);

  const std::list<std::tuple<Var, Var, double>>&
  weightSplits() const;

  const std::list<std::tuple<Var, Var, double>>&
  randomSplits() const;

  template <typename Stream>
  void
  toXML(Stream&) const;

private:
  double
  logPartSum() const;

  std::vector<std::tuple<Var, Var, double>>
  candidateParentsSplits(const Set&, const OptimalBeta&) const;

  template <typename Generator, typename SplitIt>
  bool
  chooseSplits(Generator&, const std::vector<std::tuple<Var, Var, double>>&&, const uint32_t, SplitIt, SplitIt) const;

private:
  std::list<std::tuple<Var, Var, double>> m_weightSplits;
  std::list<std::tuple<Var, Var, double>> m_randomSplits;
  std::pair<std::shared_ptr<TreeNode<Data, Var, Set>>, std::shared_ptr<TreeNode<Data, Var, Set>>> m_children;
  const Data& m_data;
  const Set m_observations;
  double m_score;
  double m_sum;
  double m_sum2;
  uint32_t m_count;
  bool m_leaf;
}; // class TreeNode

template <typename Data, typename Var, typename Set>
TreeNode<Data, Var, Set>::TreeNode(
  const Data& data,
  const Set& variables,
  const Set& observations
) : m_weightSplits(),
    m_randomSplits(),
    m_children(),
    m_data(data),
    m_observations(observations),
    m_score(0.0),
    m_sum(0.0),
    m_sum2(0.0),
    m_count(0),
    m_leaf(true)
{
  for (const auto v : variables) {
    for (const auto o : observations) {
      auto d = data(v, o);
      if (!std::isnan(d)) {
        m_sum += d;
        m_sum2 += d * d;
        ++m_count;
      }
    }
  }
  m_score = computeLogLikelihood(m_count, m_sum, m_sum2);
}

template <typename Data, typename Var, typename Set>
TreeNode<Data, Var, Set>::TreeNode(
  const std::shared_ptr<TreeNode<Data, Var, Set>>& leftChild,
  const std::shared_ptr<TreeNode<Data, Var, Set>>& rightChild
) : m_weightSplits(),
    m_randomSplits(),
    m_children(std::make_pair(leftChild, rightChild)),
    m_data(leftChild->m_data),
    m_observations(set_union(leftChild->m_observations, rightChild->m_observations)),
    m_score(0.0),
    m_sum(leftChild->m_sum + rightChild->m_sum),
    m_sum2(leftChild->m_sum2 + rightChild->m_sum2),
    m_count(leftChild->m_count + rightChild->m_count),
    m_leaf(false)
{
  m_score = computeLogLikelihood(m_count, m_sum, m_sum2);
}

template <typename Data, typename Var, typename Set>
const std::pair<std::shared_ptr<TreeNode<Data, Var, Set>>, std::shared_ptr<TreeNode<Data, Var, Set>>>&
TreeNode<Data, Var, Set>::children(
) const
{
  return m_children;
}

template <typename Data, typename Var, typename Set>
void
TreeNode<Data, Var, Set>::makeLeaf(
)
{
  m_leaf = true;
  m_children.first.reset();
  m_children.second.reset();
}

template <typename Data, typename Var, typename Set>
uint32_t
TreeNode<Data, Var, Set>::nodeCount(
  const bool includeLeaves
) const
{
  uint32_t count = 0;
  if (!m_leaf || includeLeaves) {
    ++count;
  }
  if (!m_leaf) {
    count += m_children.first->nodeCount(includeLeaves);
    count += m_children.second->nodeCount(includeLeaves);
  }
  return count;
}

template <typename Data, typename Var, typename Set>
std::list<TreeNode<Data, Var, Set>*>
TreeNode<Data, Var, Set>::nodes(
  const bool includeLeaves
)
{
  std::list<TreeNode<Data, Var, Set>*> nodes;
  if (!m_leaf || includeLeaves) {
    nodes.push_back(this);
  }
  if (!m_leaf) {
    auto firstNodes = m_children.first->nodes(includeLeaves);
    nodes.insert(nodes.end(), firstNodes.begin(), firstNodes.end());
    auto secondNodes = m_children.second->nodes(includeLeaves);
    nodes.insert(nodes.end(), secondNodes.begin(), secondNodes.end());
  }
  return nodes;
}

template <typename Data, typename Var, typename Set>
const Set&
TreeNode<Data, Var, Set>::observations(
) const
{
  return m_observations;
}

template <typename Data, typename Var, typename Set>
double
TreeNode<Data, Var, Set>::score(
) const
{
  return m_score;
}

template <typename Data, typename Var, typename Set>
double
TreeNode<Data, Var, Set>::mergeScore(
  const bool scoreBHC
) const
{
  auto mergeScore = 0.0;
  if (!m_leaf) {
    mergeScore = m_score;
    if (scoreBHC) {
      mergeScore -= (m_children.first->logPartSum() + m_children.second->logPartSum());
    }
    else {
      mergeScore -= (m_children.first->score() + m_children.second->score());
    }
  }
  return mergeScore;
}

template <typename Data, typename Var, typename Set>
double
TreeNode<Data, Var, Set>::mean(
) const
{
  return m_sum / m_count;
}

template <typename Data, typename Var, typename Set>
void
TreeNode<Data, Var, Set>::prune(
  const double scoreGain
)
{
  if (!m_leaf) {
    auto mergeScore = m_children.first->score() + m_children.second->score() -
                      m_score;
    if (mergeScore > scoreGain) {
      m_children.first->prune(scoreGain);
      m_children.second->prune(scoreGain);
    }
    else {
      LOG_MESSAGE(debug, "Pruning subtrees with merge score %g", mergeScore);
      m_children.first.reset();
      m_children.second.reset();
      m_leaf = true;
    }
  }
}

template <typename Data, typename Var, typename Set>
double
TreeNode<Data, Var, Set>::logPartSum(
) const
{
  auto lps = m_score;
  if (!m_leaf) {
    lps += log(1 + exp(m_children.first->logPartSum() + m_children.second->logPartSum() -
                       m_score));
  }
  return lps;
}

template <typename Data, typename Var, typename Set>
uint64_t
TreeNode<Data, Var, Set>::splitWeight(
  const Set& candidateParents
) const
{
  uint64_t splitCount = candidateParents.size() * m_observations.size();
  uint64_t splitWeight = splitCount * m_observations.size();
  return splitWeight;
}

template <typename Data, typename Var, typename Set>
std::vector<std::tuple<Var, Var, double>>
TreeNode<Data, Var, Set>::candidateParentsSplits(
  const Set& candidateParents,
  const OptimalBeta& ob
) const
{
  std::list<std::tuple<Var, Var, double>> splits;
  for (const auto v : candidateParents) {
    Assignment<Data, Var, Set> assmt(m_data, v, this);
    for (const auto o : m_observations) {
      auto sv = m_data(v, o);
      if (!std::isnan(sv)) {
        LOG_MESSAGE(trace, "Considering split (%s, %g)", m_data.varName(v), sv);
        auto sign = assmt.sign(sv);
        if (sign == 0) {
          LOG_MESSAGE(trace, "Split sign is zero. Skipping");
          continue;
        }
        auto beta = ob.find(assmt, sv, sign);
        if (std::isnan(beta)) {
          continue;
        }
        auto score = assmt.score(sv, sign, beta);
        LOG_MESSAGE(trace, "Parent split details for (%s, %g) : beta=%g, score=%g", m_data.varName(v), sv, beta, score);
        splits.push_back(std::make_tuple(v, o, score));
      }
    }
  }
  return std::vector<std::tuple<Var, Var, double>>(splits.begin(), splits.end());
}

template <typename Data, typename Var, typename Set>
std::list<std::tuple<Var, Var, double>>
TreeNode<Data, Var, Set>::candidateParentsSplits(
  const Set& candidateParents,
  const OptimalBeta& ob,
  const uint64_t firstWeight,
  const uint64_t maxWeight
) const
{
  std::list<std::tuple<Var, Var, double>> splits;
  auto unitWeight = m_observations.size();
  uint64_t firstSplit = (firstWeight / unitWeight) + ((firstWeight % unitWeight) ? 1 : 0);
  auto lastWeight = firstWeight + maxWeight;
  uint64_t lastSplit = (lastWeight / unitWeight) + ((lastWeight % unitWeight) ? 1 : 0) - 1;
  auto firstParent = firstSplit / m_observations.size();
  auto lastParent = lastSplit / m_observations.size();
  auto pIt = std::next(candidateParents.begin(), firstParent);
  auto prevSplits = firstSplit;
  for (auto p = firstParent; p <= lastParent; ++p, ++pIt) {
    Assignment<Data, Var, Set> assmt(m_data, *pIt, this);
    auto firstObservation = prevSplits % m_observations.size();
    auto numObservations = std::min(m_observations.size() - firstObservation, lastSplit + 1 - prevSplits);
    auto oIt = std::next(m_observations.begin(), firstObservation);
    for (auto o = firstObservation; o < firstObservation + numObservations; ++o, ++oIt) {
      auto sv = m_data(*pIt, *oIt);
      if (!std::isnan(sv)) {
        LOG_MESSAGE(trace, "Considering split (%s, %g)", m_data.varName(*pIt), sv);
        auto sign = assmt.sign(sv);
        if (sign == 0) {
          LOG_MESSAGE(trace, "Split sign is zero. Skipping");
          continue;
        }
        auto beta = ob.find(assmt, sv, sign);
        if (std::isnan(beta)) {
          continue;
        }
        auto score = assmt.score(sv, sign, beta);
        LOG_MESSAGE(trace, "Parent split details for (%s, %g) : beta=%g, score=%g", m_data.varName(*pIt), sv, beta, score);
        splits.push_back(std::make_tuple(*pIt, *oIt, score));
      }
    }
    prevSplits += numObservations;
  }
  return splits;
}

template <typename Data, typename Var, typename Set>
template <typename Generator, typename SplitIt>
bool
TreeNode<Data, Var, Set>::chooseSplits(
  Generator& generator,
  const std::vector<std::tuple<Var, Var, double>>&& candidateSplits,
  const uint32_t numSplits,
  SplitIt weightIt,
  SplitIt randomIt
) const
{
  if (candidateSplits.empty()) {
    // Move generator state forward to simulate picking
    // 2 * numSplits splits
    ::advance(generator, 2 * numSplits);
    LOG_MESSAGE(debug, "No candidate splits found");
    return false;
  }
  LOG_MESSAGE(debug, "Number of candidate splits found: %u", candidateSplits.size());
  auto maxScoreSplit = std::max_element(candidateSplits.cbegin(), candidateSplits.cend(),
                                        [] (const std::tuple<Var, Var, double>& a,
                                            const std::tuple<Var, Var, double>& b)
                                           { return std::isless(std::get<2>(a), std::get<2>(b)); });
  auto maxScore = std::get<2>(*maxScoreSplit);
  std::vector<double> weights(candidateSplits.size());
  std::transform(candidateSplits.cbegin(), candidateSplits.cend(), weights.begin(),
                 [&maxScore] (const std::tuple<Var, Var, double>& s)
                             { return exp(std::get<2>(s) - maxScore); });
  discrete_distribution_safe<uint64_t> splitWeight(weights.cbegin(), weights.cend());
  trng::uniform_int_dist splitRand(0, candidateSplits.size());
  for (auto i = 0u; i < numSplits; ++i, ++weightIt, ++randomIt) {
    auto w = splitWeight(generator);
    *weightIt = candidateSplits[w];
    LOG_MESSAGE(debug, "Chosen parent split using weights: (%s, %g)",
                       m_data.varName(std::get<0>(*weightIt)),
                       m_data(std::get<0>(*weightIt), std::get<1>(*weightIt)));
    auto r = splitRand(generator);
    *randomIt = candidateSplits[r];
    LOG_MESSAGE(debug, "Chosen parent split at random: (%s, %g)",
                       m_data.varName(std::get<0>(*randomIt)),
                       m_data(std::get<0>(*randomIt), std::get<1>(*randomIt)));
  }
  return true;
}

template <typename Data, typename Var, typename Set>
template <typename Generator>
bool
TreeNode<Data, Var, Set>::learnParentsSplits(
  Generator& generator,
  const Set& candidateParents,
  const OptimalBeta& ob,
  const uint32_t numSplits
)
{
  auto splits = this->candidateParentsSplits(candidateParents, ob);
  if (!splits.empty()) {
    m_weightSplits.resize(numSplits);
    m_randomSplits.resize(numSplits);
  }
  return this->chooseSplits(generator, std::move(splits), numSplits, m_weightSplits.begin(), m_randomSplits.begin());
}

template <typename Data, typename Var, typename Set>
template <typename Generator, typename SplitIt>
bool
TreeNode<Data, Var, Set>::learnParentsSplits(
  Generator& generator,
  const Set& candidateParents,
  const OptimalBeta& ob,
  const uint32_t numSplits,
  SplitIt weightIt,
  SplitIt randomIt
) const
{
  auto splits = this->candidateParentsSplits(candidateParents, ob);
  return this->chooseSplits(generator, std::move(splits), numSplits, weightIt, randomIt);
}

template <typename Data, typename Var, typename Set>
template <typename SplitIt>
void
TreeNode<Data, Var, Set>::setWeightSplits(
  const SplitIt& first,
  const SplitIt& last
)
{
  m_weightSplits.clear();
  std::copy(first, last, std::back_inserter(m_weightSplits));
}

template <typename Data, typename Var, typename Set>
template <typename SplitIt>
void
TreeNode<Data, Var, Set>::setRandomSplits(
  const SplitIt& first,
  const SplitIt& last
)
{
  m_randomSplits.clear();
  std::copy(first, last, std::back_inserter(m_randomSplits));
}

template <typename Data, typename Var, typename Set>
const std::list<std::tuple<Var, Var, double>>&
TreeNode<Data, Var, Set>::weightSplits(
) const
{
  return m_weightSplits;
}

template <typename Data, typename Var, typename Set>
const std::list<std::tuple<Var, Var, double>>&
TreeNode<Data, Var, Set>::randomSplits(
) const
{
  return m_randomSplits;
}

template <typename Data, typename Var, typename Set>
template <typename Stream>
void
TreeNode<Data, Var, Set>::toXML(
  Stream& stream
) const
{
  if (m_leaf) {
    stream << "<TreeNode numChildren=\"0\" score=\"" << m_score <<
              "\" condSet=\"" << m_observations << "\">";
  }
  else {
    stream << "<TreeNode numChildren=\"2\" score=\"" << m_score <<
              "\" condSet=\"" << m_observations << "\">" << std::endl;
    for (const auto& split : m_weightSplits) {
      stream << "<Regulator name=\"" << m_data.varName(std::get<0>(split)) <<
                "\" entropy=\"" << std::get<2>(split) <<
                "\" splitValue=\"" << m_data(std::get<0>(split), std::get<1>(split)) <<
                "\">";
      stream << "</Regulator>" << std::endl;
    }
    for (const auto& split : m_randomSplits) {
      stream << "<RandomRegulator name=\"" << m_data.varName(std::get<0>(split)) <<
                "\" entropy=\"" << std::get<2>(split) <<
                "\" splitValue=\"" << m_data(std::get<0>(split), std::get<1>(split)) <<
                "\">";
      stream << "</RandomRegulator>" << std::endl;
    }
    m_children.first->toXML(stream);
    m_children.second->toXML(stream);
  }
  stream << "</TreeNode>" << std::endl;
}

#endif // DETAIL_TREENODE_HPP_
