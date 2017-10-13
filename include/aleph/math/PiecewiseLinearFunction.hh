#ifndef ALEPH_MATH_PIECEWISE_LINEAR_FUNCTION_HH__
#define ALEPH_MATH_PIECEWISE_LINEAR_FUNCTION_HH__

#include <iterator>
#include <map>
#include <set>
#include <stdexcept>

namespace aleph
{

namespace math
{

template <class D, class I = D> class PiecewiseLinearFunction
{
public:
  using Domain = D;
  using Image  = I;

  /**
    Creates a new piecewise linear function from a range of values. The values
    must consist of pairs that are convertible to double. Duplicates are not
    permitted and will result in an exception being thrown.

    @param begin Input iterator to begin of range
    @param end   Input iterator to end of range
  */

  template <class InputIterator> PiecewiseLinearFunction( InputIterator begin,
                                                          InputIterator end )
  {
    for( InputIterator it = begin; it != end; ++it )
    {
      auto x       = it->first;
      auto y       = it->second;
      bool success = false;

      std::tie( std::ignore, success )
          = _data.insert( std::make_pair(x,y) );

      if( !success )
        throw std::runtime_error( "Duplicate value pairs not permitted for piecewise linear functions" );
    }

    // Inserts new points whenever the function intersects the x-axis in
    // order to ensure that the integration procedure works.
    if( !_data.empty() )
      this->insertIntersectionPoints();
  }

private:

  /**
    Checks the segments of the piecewise linear function for intersections with
    the x-axis. In these cases, a new point needs to be inserted into the list
    of function values. Else, the integration procedure will not work because a
    segment might be considered to have an area of zero. The following example
    illustrates this problem:

    \verbatim
     /\
    /  \
    ---------
        \  /
         \/
    \endverbatim

    If the intersection in the middle is not being counted as a point on the PL
    function, the area of the corresponding segment will be zero. This is not
    the desired and expected behaviour.
  */

  void insertIntersectionPoints()
  {
    auto current = _data.begin();
    auto next    = std::next( current );

    std::set<Domain> intersections;

    for( ; next != _data.end(); )
    {
      auto x0 = current->first;
      auto y0 = current->second;

      auto x1 = next->first;
      auto y1 = next->second;

      // We do not need to check the other cases. If either one of the values is
      // zero, we already have an intersection.
      if( y0 * y1 < Image() )
      {
        auto m = ( y1 - y0 ) / ( x1 - x0 );
        auto c = y0 - m * x0;
        auto x = - c / m;

        intersections.insert( x );
      }

      current = next;
      next    = std::next( current );
    }

    for( auto&& x : intersections )
      _data.insert( std::make_pair( x, Image() ) );
  }

  /** Maps values in the domain to values in the image */
  std::map<Domain, Image> _data;
};

} // namespace math

} // namespace aleph

#endif
