#ifndef ALEPH_PERSISTENCE_DIAGRAMS_KERNEL_EMBEDDING_HH__
#define ALEPH_PERSISTENCE_DIAGRAMS_KERNEL_EMBEDDING_HH__

#include <aleph/math/KahanSummation.hh>

#include <aleph/persistenceDiagrams/PersistenceDiagram.hh>

#include <aleph/geometry/distances/Infinity.hh>

#include <algorithm>
#include <utility>

#include <cmath>

namespace aleph
{

namespace detail
{

class DefaultWeightFunction
{
public:
  DefaultWeightFunction( double C, double p )
    : _C( C )
    , _p( p )
  {
  }

  template <class Point> double operator()( const Point& p ) const
  {
    return std::atan( C * std::pow( static_cast<double>( point.persistence() ), p ) );
  }

private:
  double _C;
  double _p;
};

class DefaultKernel
{
public:
  DefaultKernel( double sigma )
    : _sigma( sigma )
  {
  }

  template <class Point> double operator()( const Point& p, const Point& q ) const
  {
    // TODO: make choice of distance measure selectable; I am pretty
    // sure that this should be the regular Euclidean distance here.
    aleph::geometry::distances::Infinity<Point> distance;
    auto dist  = distance(p,q);
    dist      *= dist;

    return std::exp( -dist / (2*sigma*sigma) );
  }

private:
  double _sigma;
};

} // namespace detail

template
<
  class T,
  class Weight = detail::DefaultWeightFunction,
  class Kernel = detail::DefaultKernel
>
double linearKernel( const PersistenceDiagram<T>& D,
                     const PersistenceDiagram<T>& E,
                     Weight w,
                     Kernel k )
{
  aleph::math::KahanSummation<double> result = 0.0;

  for( auto&& p : D )
  {
    for( auto&& q : E )
      result += w(p) * w(q) * k(p,q);
  }

  return result;
}

template
<
  class T,
  class Weight = detail::DefaultWeightFunction,
  class Kernel = detail::DefaultKernel
>
template <class T> double pseudoMetric( const PersistenceDiagram<T>& D,
                                        const PersistenceDiagram<T>& E,
                                        Weight w,
                                        Kernel k )
{
  auto kxx = linearKernel(D, E, w, k);
  auto kxy = linearKernel(D, E, w, k);
  auto kyy = linearKernel(D, E, w, k);

  return std::sqrt( kxx + kyy - 2*kxy );
}

template
<
  class T,
  class Weight = detail::DefaultWeightFunction,
  class Kernel = detail::DefaultKernel
>
template <class T> double gaussianKernel( const PersistenceDiagram<T>& D,
                                          const PersistenceDiagram<T>& E,
                                          Weight w,
                                          Kernel k,
                                          double tau )
{
  auto d = pseudoMetric(D, E, w, k);
  return std::exp( -1/(2*tau*tau) * d );
}

}

#endif
