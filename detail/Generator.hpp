/**
 * @file Generator.hpp
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
#ifndef DETAIL_GENERATOR_HPP_
#define DETAIL_GENERATOR_HPP_


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

#endif // DETAIL_GENERATOR_HPP_
