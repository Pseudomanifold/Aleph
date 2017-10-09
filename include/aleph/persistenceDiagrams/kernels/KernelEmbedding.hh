#ifndef ALEPH_PERSISTENCE_DIAGRAMS_KERNEL_EMBEDDING_HH__
#define ALEPH_PERSISTENCE_DIAGRAMS_KERNEL_EMBEDDING_HH__

#include <aleph/math/KahanSummation.hh>

#include <aleph/persistenceDiagrams/PersistenceDiagram.hh>

#include <aleph/geometry/distances/Euclidean.hh>

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
    aleph::geometry::distances::Euclidean<double> distance;
    auto dist  = distance(p,q);
    return std::exp( -dist / (2*sigma*sigma) );
  }

private:
  double _sigma;
};

} // namespace detail

/**
  Calculates the linear version of the persistence-weighted Gaussian
  kernel between two persistence diagrams.

  @param D First persistence diagram
  @param E Second persistence diagram
  @param w Weight function; by default, a function based on `atan` is being used
  @param k Kernel function; by default, a *Gaussian kernel* is used. Notice that
           the parameters of this kernel need to be adjusted.

  @tparam T      Data type of persistence diagram (inferred)
  @tparam Weight Functor for calculating weights of persistence points
  @tparam Kernel Functor for calculating kernel values of persistence points

  @returns Kernel value
*/

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

/**
  Calculate the pseudo-metric based on the persistence-weighted Gaussian
  kernel. A linear kernel is used to obtain a value for the metric. This
  follows the approach in the original paper.

  @see http://proceedings.mlr.press/v48/kusano16.pdf (original paper by Kusano et al.)

  @param D First persistence diagram
  @param E Second persistence diagram
  @param w Weight function; by default, a function based on `atan` is being used
  @param k Kernel function; by default, a *Gaussian kernel* is used. Notice that
           the parameters of this kernel need to be adjusted.

  @tparam T      Data type of persistence diagram (inferred)
  @tparam Weight Functor for calculating weights of persistence points
  @tparam Kernel Functor for calculating kernel values of persistence points

  @returns Pseudo-metric value
*/

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

/**
  Calculate the Gaussian kernel value based on the persistence-weighted
  Gaussian kernel, using a smoothing parameter \p tau. This function is
  using the *pseudo-metric function*, which in turn employs a *linear*
  kernel.

  @see http://proceedings.mlr.press/v48/kusano16.pdf (original paper by Kusano et al.)

  @param D   First persistence diagram
  @param E   Second persistence diagram
  @param w   Weight function; by default, a function based on `atan` is being used
  @param k   Kernel function; by default, a *Gaussian kernel* is used. Notice that
             the parameters of this kernel need to be adjusted.
  @param tau Smoothing parameter

  @tparam T      Data type of persistence diagram (inferred)
  @tparam Weight Functor for calculating weights of persistence points
  @tparam Kernel Functor for calculating kernel values of persistence points

  @returns Gaussian kernel value
*/

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

} // namespace aleph

#endif
