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
  Scalar multiplication of a vector.

  @param lambda Scalar
  @param vector Vector

  @returns Vector in which every element has been multiplied by the scalar
*/

template <class U, class V> std::vector<U> operator*( V lambda, const std::vector<U>& vector )
{
  std::vector<U> result( vector.size() );

  std::transform( vector.begin(), vector.end(), result.begin(),
    [&lambda] ( U x )
    {
      return static_cast<U>(lambda) * x;
    }
  );

  return result;
}

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

/**
  Checks whether all elements of two sequences are close to each other
  with some pre-defined tolerance. This function *closely* follows the
  `numpy.allclose` function.
*/

template <class InputIterator1, class InputIterator2> bool allclose(
  InputIterator1 begin1, InputIterator1 end1,
  InputIterator2 begin2, InputIterator2 end2,
  double rtol = 1e-05,
  double atol = 1e-08 )
{
  auto n = std::distance( begin1, end1 );
  auto m = std::distance( begin2, end2 );

  if( n != m )
    return false;

  auto it1 = begin1;
  auto it2 = begin2;

  for( ; it1 != end1 && it2 != end2; ++it1, ++it2 )
  {
    if( double( std::abs( *it1 -  *it2 ) ) > ( atol + rtol * std::abs( double( *it2 ) ) ) )
      return false;
  }

  return true;
}

} // namespace utilities

} // namespace aleph

#endif
