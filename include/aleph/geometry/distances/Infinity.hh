#ifndef ALEPH_GEOMETRY_DISTANCES_INFINITY_HH__
#define ALEPH_GEOMETRY_DISTANCES_INFINITY_HH__

#include <algorithm>

namespace aleph
{

namespace geometry
{

namespace distances
{

/**
  Basic functor for calculating the infinity distance between two
  points. The point class needs to provide coordinate access. The
  functor is mainly used for persistence diagram distances.
*/

template <class DataType> class InfinityDistance
{
public:
  template <class Point> DataType operator()( const Point& p, const Point& q ) const
  {
    // This ensures that the values stay positive, regardless of the
    // underlying data type.
    auto dx = p.x() >= q.x() ? p.x() - q.x() : q.x() - p.x();
    auto dy = p.y() >= q.y() ? p.y() - q.y() : q.y() - p.y();

    return std::max( dx, dy );
  }
};

} // namespace distances

} // namespace geometry

} // namespace aleph

#endif
