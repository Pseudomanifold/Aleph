#ifndef ALEPH_UTILITIES_CONTAINER_OPERATORS_HH__
#define ALEPH_UTILITIES_CONTAINER_OPERATORS_HH__

/**
  @file  ContainerOperators.hh
  @brief Defines (arithmetic) operators to use on standard containers

  In some cases, it makes sense to treat an STL container as
  a mathematical vector. These objects then should support a
  set of arithmetic operations, such as addition. As the STL
  does not provide these things, we have to do them manually
  in this file.

  Note that the operations are kept in *aleph::utilities* in
  order not to pollute the global namespace.
*/

#include <algorithm>
#include <iterator>
#include <vector>

#include <cassert>

namespace aleph
{

namespace utilities
{

/**
  Element-wise addition of two vectors of the same size.

  @param lhs First vector
  @param rhs Second vector

  @returns Addition of the two vectors
*/

template <class T> std::vector<T> operator+( const std::vector<T>& lhs, const std::vector<T>& rhs )
{
  assert( lhs.size() == rhs.size() );

  std::vector<T> result;
  result.reserve( lhs.size() );

  std::transform( lhs.begin(), lhs.end(), rhs.begin(), std::back_inserter( result ), std::plus<T>() );
  return result;
}

/**
  Element-wise subtraction of two vectors of the same size.

  @param lhs First vector
  @param rhs Second vector

  @returns Addition of the two vectors
*/

template <class T> std::vector<T> operator-( const std::vector<T>& lhs, const std::vector<T>& rhs )
{
  assert( lhs.size() == rhs.size() );

  std::vector<T> result;
  result.reserve( lhs.size() );

  std::transform( lhs.begin(), lhs.end(), rhs.begin(), std::back_inserter( result ), std::minus<T>() );
  return result;
}

} // namespace utilities

} // namespace aleph

#endif
