#ifndef ALEPH_DISTANCES_INFINITY_HH__
#define ALEPH_DISTANCES_INFINITY_HH__

#include "PersistenceDiagram.hh"

#include <algorithm>

namespace aleph
{

namespace distances
{

template <class DataType> class InfinityDistance
{
public:
  using PersistenceDiagram = PersistenceDiagram<DataType>;
  using Point              = typename PersistenceDiagram::Point;

  DataType operator()( const Point& p, const Point& q ) const
  {
    // This ensures that the values stay positive, regardless of the
    // underlying data type.
    auto dx = p.x() >= q.x() ? p.x() - q.x() : q.x() - p.x();
    auto dy = p.y() >= q.y() ? p.y() - q.y() : q.y() - p.y();

    return std::max( dx, dy );
  }
};

}

}

#endif
