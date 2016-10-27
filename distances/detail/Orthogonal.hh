#ifndef ALEPH_DISTANCES_DETAIL_ORTHOGONAL_HH__
#define ALEPH_DISTANCES_DETAIL_ORTHOGONAL_HH__

namespace aleph
{

namespace distances
{

namespace detail
{

template <class Distance, class Point> auto orthogonalDistance( const Point& p ) -> decltype( Distance().operator()( p,p ) )
{
  auto&& x = p.x();
  auto&& y = p.y();

  // Here's to hoping that the division by two is not being truncated by the
  // underlying data type. If it is, there's nothing we can do...
  auto&& u = ( x + y ) / 2;

  Point q( u,u );
  return Distance().operator()( p,q );
}

}

}

}

#endif
