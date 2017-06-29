#ifndef ALEPH_GEOMETRY_DISTANCES_MANHATTAN_HH__
#define ALEPH_GEOMETRY_DISTANCES_MANHATTAN_HH__

#include <cstddef>
#include <cmath>

#include <iterator>
#include <string>

namespace aleph
{

namespace distances
{

/**
  Manhattan distance ($L_1$ distance) functor The functor is templated in order
  to handle arbitrary ranges. Its interface is general, but the \c accum_dist
  function is specifically meant to be used with FLANN.

  @see FLANN project page (http://www.cs.ubc.ca/research/flann)
  @see FLANN repository   (https://github.com/mariusmuja/flann)
*/

template <class T> class Manhattan
{
public:

  // Flag telling FLANN that this functor can be used for calculating distances
  // within kd trees.
  using is_kdtree_distance = bool;

  // Required for FLANN usage
  using ElementType = T;
  using ResultType  = T;

  /**
    Given two ranges of double values, which are assumed to represent two
    vectors, calculates the distance between them.

    @param a             Iterator describing first vector
    @param b             Iterator describing second vector
    @param size          Size of vectors a and b

    @param worstDistance If set to a value greater than zero, calculations will
    stop once the value has been reached. Else, the value of this variable is
    ignored. This variable is provided in order to remain compatible with the
    interface of FLANN.

    @returns Manhattan distance between the two input vectors.
  */

  template <typename Iterator1, typename Iterator2>
  ResultType operator()( Iterator1 a,
                         Iterator2 b,
                         std::size_t size,
                         ElementType worstDistance = -1.0 ) const
  {
    // Fixes warnings about unused parameters. This parameter is
    // provided for compatibility reasons with FLANN only.
    (void) worstDistance;

    ElementType result = 0.0;

    ElementType diff0  = 0.0;
    ElementType diff1  = 0.0;
    ElementType diff2  = 0.0;
    ElementType diff3  = 0.0;

    using DifferenceType1
      = typename std::iterator_traits<Iterator1>::difference_type;

    Iterator1 last      = std::next( a, DifferenceType1( size ) );
    Iterator1 lastGroup = std::prev( last, 3 );

    while( a < lastGroup )
    {
      diff0 = std::abs( ElementType( a[0] - b[0] ) );
      diff1 = std::abs( ElementType( a[1] - b[1] ) );
      diff2 = std::abs( ElementType( a[2] - b[2] ) );
      diff3 = std::abs( ElementType( a[3] - b[3] ) );

      result +=   diff0
                + diff1
                + diff2
                + diff3;

      a += 4;
      b += 4;

      if( worstDistance > 0 && result > worstDistance )
        return result;
    }

    while( a < last )
    {
      diff0  = std::abs( ElementType( *a++ - *b++ ) );
      result = result + diff0;
    }

    return result;
  }

  /**
    Partial distance calculation, used by FLANN for fast kd-tree calculations.
    This function exploits that the Manhattan distance can be evaluated
    component-wise.

    @param a First component
    @param b Second component

    @returns Partial distance between those two components
  */

  template <typename U, typename V>
  ResultType accum_dist( const U& a,
                         const V& b,
                         int __attribute__((unused)) ) const
  {
    return std::abs( a - b );
  }

  /** @returns Name of functor */
  static std::string name()
  {
    return "Manhattan distance";
  }
};

} // namespace distances

} // namespace aleph

#endif
