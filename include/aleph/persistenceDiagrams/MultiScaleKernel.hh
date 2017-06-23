#ifndef ALEPH_MULTI_SCALE_KERNEL_HH__
#define ALEPH_MULTI_SCALE_KERNEL_HH__

#include <aleph/math/KahanSummation.hh>

#include <aleph/persistenceDiagrams/PersistenceDiagram.hh>

#include <algorithm>
#include <utility>

#include <cmath>

namespace aleph
{

namespace detail
{

/**
  Calculates the squared Euclidean distance between points in
  a persistence diagram. If specified, the second point, $q$,
  is mirrored at the diagonal.
*/

template <class Point> double squaredEuclideanDistance( const Point& p,
                                                        const Point& q,
                                                        bool mirror = false )
{
  auto x0 = p.x();
  auto y0 = p.y();
  auto x1 = q.x();
  auto y1 = q.y();

  if( mirror )
    std::swap( x1, y1 );

  // The check ensures that dx and dy are always positive. Thus, even if
  // the underlying data type is an unsigned type, both expressions will
  // evaluate properly.
  auto dx = x0 > x1 ? x0-x1 : x1-x0;
  auto dy = y0 > y1 ? y0-y1 : y1-y0;

  return static_cast<double>( dx*dx + dy*dy );
}

} // namespace detail

template <class T> double multiScaleKernel( const PersistenceDiagram<T>& D1,
                                            const PersistenceDiagram<T>& D2,
                                            double sigma )
{
  aleph::math::KahanSummation<double> sum = 0.0;

  for( auto&& p : D1 )
  {
    for( auto&& q : D2 )
    {
      auto d1 = detail::squaredEuclideanDistance( p, q );
      auto d2 = detail::squaredEuclideanDistance( p, q, true );

      sum += std::exp( -d1 / ( 8.0*M_PI ) );
      sum -= std::exp( -d2 / ( 8.0*M_PI ) );
    }
  }

  return 1.0 / ( 8.0*M_PI*sigma ) * sum;
}

template <class T> double multiScaleFeatureMap( const PersistenceDiagram<T>& D,
                                                double sigma )
{
  aleph::math::KahanSummation<double> result = 0.0;

  for( auto&& p : D )
  {
    for( auto&& q : D )
    {
      auto d1 = detail::squaredEuclideanDistance( p, q );
      auto d2 = detail::squaredEuclideanDistance( p, q, true );

      result += std::exp( -d1 / (4.0*M_PI) );
      result -= std::exp( -d2 / (4.0*M_PI) );
    }
  }

  return 1.0 / ( 4.0*M_PI*sigma ) * result;
}

template <class T> double multiScalePseudoMetric( const PersistenceDiagram<T>& D1,
                                                  const PersistenceDiagram<T>& D2,
                                                  double sigma )
{
  auto kxx = multiScaleKernel( D1, D1, sigma );
  auto kxy = multiScaleKernel( D1, D2, sigma );
  auto kyy = multiScaleKernel( D2, D2, sigma );

  return std::sqrt( kxx + kyy - 2*kxy );
}

}

#endif
