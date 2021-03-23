/**
 * @file Random.hpp
 * @brief Implementation of some PRNG functions.
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
#ifndef DETAIL_RANDOM_HPP_
#define DETAIL_RANDOM_HPP_

#include <trng/discrete_dist.hpp>

#include <algorithm>


/**
 * @brief Safe version of trng::fast_discrete_dist
 *        which handles infinite weights.
 *
 * @tparam IntType Type of integer for picking the index.
 */
template <typename IntType = int>
class discrete_distribution_safe {
public:
 /**
  * @brief Constructor for the distribution.
  *
  * @tparam InputIt Type of input iterator with weights.
  * @param first Iterator to the first weight.
  * @param last Iterator to the last weight.
  */
  template <typename InputIt>
  discrete_distribution_safe(
    const InputIt first,
    const InputIt last
  ) : m_distrib(first, last),
      m_infiniteIndex(std::numeric_limits<IntType>::max()),
      m_infiniteWeight(false)
  {
    static auto findInf = [] (const double w) { return std::isinf(w); };
    auto it = std::find_if(first, last, findInf);
    if (it != last) {
      m_infiniteIndex = std::distance(first, it);
      m_infiniteWeight = true;
    }
  }

  /**
   * @brief Returns a random index based on the given weights if all
   *        the weights are finite. Otherwise, returns the index of
   *        the first infinite weight. Always moves PRNG state forward.
   *
   * @tparam GenType Type of the pseudo random number generator.
   * @param generator Instance of the PRNG.
   */
  template <typename Generator>
  IntType
  operator()(Generator& generator)
  {
    if (!m_infiniteWeight) {
      return static_cast<IntType>(m_distrib(generator));
    }
    else {
      generator.discard(1);
      return m_infiniteIndex;
    }
  }

private:
  trng::discrete_dist m_distrib;
  IntType m_infiniteIndex;
  bool m_infiniteWeight;
}; // class discrete_distribution_safe


template <typename Generator>
class HasJump {
private:
    typedef uint8_t YesType;
    typedef uint16_t NoType;

private:
    template <typename G>
    static
    YesType&
    test(decltype(&G::jump));

    template <typename G>
    static
    NoType&
    test(...);

public:
    enum {
      value = (sizeof(test<Generator>(0)) == sizeof(YesType))
    };
}; // class HasJump<Generator>

template <typename Generator>
typename std::enable_if<HasJump<Generator>::value, void>::type
advance(
  Generator& generator,
  const uint64_t count
)
{
  generator.jump(count);
}

template <typename Generator>
typename std::enable_if<!HasJump<Generator>::value, void>::type
advance(
  Generator& generator,
  const uint64_t count
)
{
  generator.discard(count);
}

#endif // DETAIL_RANDOM_HPP_
