#ifndef ALEPH_GEOMETRY_DISTANCES_EUCLIDEAN_HH__
#define ALEPH_GEOMETRY_DISTANCES_EUCLIDEAN_HH__

#include <aleph/geometry/distances/Traits.hh>

#include <cmath>
#include <cstddef>

#include <iterator>
#include <string>

namespace aleph
{

namespace geometry
{

namespace distances
{

/**
  Euclidean distance ($L_2$ distance) functor. The functor is templated in
  order to handle arbitrary ranges. Its interface is somewhat general, but
  the \c accum_dist function is specifically meant to be used with FLANN.

  @see FLANN project page (http://www.cs.ubc.ca/research/flann)
  @see FLANN repository   (https://github.com/mariusmuja/flann)
*/

template <class T> class Euclidean
{
public:

  // Flag telling FLANN that this functor can be used for calculating distances
  // within kd trees.
  typedef bool is_kdtree_distance;

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
  */

  template <typename Iterator1, typename Iterator2>
  ResultType operator()( Iterator1 a,
                         Iterator2 b,
                         std::size_t size,
                         ElementType worstDistance = -1.0 ) const
  {
    // Fix compiler warnings about unused parameters. This is provided to be
    // compatible with FLANN.
    (void) worstDistance;

    ResultType result = 0.0;

    ResultType diff0  = 0.0;
    ResultType diff1  = 0.0;
    ResultType diff2  = 0.0;
    ResultType diff3  = 0.0;

    Iterator1 last      = a;
    Iterator1 lastGroup = last;

    using DifferenceType = typename std::iterator_traits<Iterator1>::difference_type;

    std::advance( last,      DifferenceType( size ) );
    std::advance( lastGroup, DifferenceType( -3 ) );

    while( a < lastGroup )
    {
      diff0 = ElementType( a[0] - b[0] );
      diff1 = ElementType( a[1] - b[1] );
      diff2 = ElementType( a[2] - b[2] );
      diff3 = ElementType( a[3] - b[3] );

      result +=   diff0 * diff0
                + diff1 * diff1
                + diff2 * diff2
                + diff3 * diff3;

      a += 4;
      b += 4;

      if( worstDistance > 0 && result > worstDistance )
        return result;
    }

    while( a < last )
    {
      diff0  = ElementType( *a++ - *b++ );
      result = result + diff0 * diff0;
    }

    return result;
  }

  /**
    Partial distance calculation, used by FLANN for fast kd-tree calculations.
    This function exploits that the Euclidean distance can be evaluated
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
    return (a-b) * (a-b);
  }

  /** @returns Name of functor */
  static std::string name()
  {
    return "Euclidean distance";
  }
};

template <class T> struct Traits< Euclidean<T> >
{
  using ResultType  = typename Euclidean<T>::ResultType;
  using ElementType = typename Euclidean<T>::ElementType;

  ResultType from( ElementType x ) const noexcept
  {
    return ResultType( std::sqrt( x ) );
  }

  ResultType to( ElementType x ) const noexcept
  {
    return ResultType( x*x );
  }
};

} // namespace geometry

} // namespace distances

} // namespace aleph

#endif
