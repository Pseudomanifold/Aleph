#ifndef ALEPH_CONTAINERS_DIMENSIONALITY_ESTIMATORS_HH__
#define ALEPH_CONTAINERS_DIMENSIONALITY_ESTIMATORS_HH__

#include <vector>

#include <aleph/math/KahanSummation.hh>

namespace aleph
{

namespace containers
{

/**
  Estimates local intrinsic dimensionality of a container using its
  nearest neighbours. The underlying assumption of the estimator is
  that points are locally uniformly distributed uniformly. Use this
  estimator with care when analysing unknown data.

  :param container: Container to use for dimensionality estimation
  :param k:         Number of nearest neighbours
  :param distance:  Distance measure

  :returns: Vector of local intrinsic dimensionality estimates that
            are reported *without* rounding.
*/

template <
  class Distance,
  class Container,
  class Wrapper
> std::vector<double> estimateLocalDimensionalityNearestNeighbours( const Container& container,
                                                                    unsigned k,
                                                                    Distance /* distance */ = Distance() )
{
  using IndexType = typename Wrapper::IndexType;

  std::vector< std::vector<IndexType> > indices;
  std::vector< std::vector<double> > distances;

  Wrapper nnWrapper( container );
  nnWrapper.neighbourSearch( k+1, indices, distances );

  auto n = container.size();

  std::vector<double> estimates;
  estimates.reserve( n );

  for( decltype(n) i = 0; i < n; i++ )
  {
    auto&& nnDistances = distances.at(i);
    auto r1            = aleph::math::accumulate_kahan( nnDistances.begin(), nnDistances.begin() + k,     0.0 ) / static_cast<double>(k  );
    auto r2            = aleph::math::accumulate_kahan( nnDistances.begin(), nnDistances.begin() + k + 1, 0.0 ) / static_cast<double>(k+1);

    estimates.push_back( r1 / ( (r2-r1)*k ) );
  }

  return estimates;

}

} // namespace containers

} // namespace aleph

#endif
