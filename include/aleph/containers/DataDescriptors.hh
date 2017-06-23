#ifndef ALEPH_CONTAINERS_DATA_DESCRIPTORS_HH__
#define ALEPH_CONTAINERS_DATA_DESCRIPTORS_HH__

#include <cmath>

#include <algorithm>
#include <vector>

#include <aleph/distances/Euclidean.hh>
#include <aleph/distances/Traits.hh>

#include <aleph/geometry/BruteForce.hh>
#include <aleph/geometry/NearestNeighbours.hh>

#include <aleph/math/KahanSummation.hh>

namespace aleph
{

/**
  Order-based eccentricity calculations. Eccentricity assigns low values
  to central points in a point cloud without having to define the actual
  centre of it. The order parameter can be used to decrease how much the
  small distances influence the result.
*/

template <class Distance, class Container> std::vector<double> eccentricities( const Container& container,
                                                                               unsigned order = 1 )
{
  auto n = container.size();
  auto d = container.dimension();

  // Ensures that subsequent algorithms always operate on a non-empty
  // range of iterators.
  if( n == 0 )
    return {};

  std::vector<double> eccentricities;
  eccentricities.reserve( n );

  Distance dist;
  distances::Traits<Distance> traits;

  for( decltype(n) i = 0; i < n; i++ )
  {
    auto p = container[i];

    std::vector<double> distances;
    distances.reserve( n );

    double eccentricity = 0.0;

    for( decltype(n) j = 0; j < n; j++ )
    {
      if( i != j )
      {
        auto q        = container[j];
        auto distance = dist( p.begin(),
                              q.begin(),
                              d );

        distance = traits.from( distance );

        if( order > 0 )
          distance = std::pow( distance, decltype(distance)( order ) ) / decltype(distance)(n);

        distances.push_back( distance );
      }
    }

    if( order > 0 )
    {
      eccentricity = aleph::math::accumulate_kahan( distances.begin(), distances.end(), 0.0 );
      eccentricity = std::pow( eccentricity, 1.0 / order );
    }
    else
      eccentricity = *std::max_element( distances.begin(), distances.end() );

    eccentricities.push_back( eccentricity );
  }

  return eccentricities;
}

/**
  Density estimation using a truncated Gaussian estimator. Points whose
  Euclidean distance is smaller than the bandwidth will not be used for
  estimating the density.
*/

template <class Container> std::vector<double> estimateDensityTruncatedGaussian( const Container& container, double bandwidth )
{
  auto n = container.size();
  auto d = container.dimension();

  const auto bandwidthSquare = bandwidth * bandwidth;

  std::vector<double> densities;
  densities.reserve( n );

  distances::Euclidean<double> distanceFunctor;

  for( decltype(n) i = 0; i < n; i++ )
  {
    auto p                                      = container[i];
    aleph::math::KahanSummation<double> density = 0.0;

    for( decltype(n) j = 0; j < n; j++ )
    {
      auto q        = container[j];
      auto distance = distanceFunctor( p.begin(),
                                       q.begin(),
                                       d );
      if( distance <= bandwidthSquare )
        density += std::exp( -1.0 * distance / ( 2.0 * bandwidth ) );
    }

    densities.push_back( density / static_cast<double>(n) );
  }

  return densities;
}

/**
  Density estimator using the distance to a measure density estimator as
  introduced by Chazal et al. in:

      Persistence-based clustering in Riemannian manifolds

  This density estimator is capable of using different distance
  functors.
*/

template <
  class Distance,
  class Container,
  class Wrapper = geometry::BruteForce<Container, Distance>
> std::vector<double> estimateDensityDistanceToMeasure( const Container& container,
                                                        unsigned k,
                                                        Distance /* distance */ = Distance() )
{
  auto n = container.size();

  std::vector<double> densities;
  densities.reserve( n );

  using IndexType = typename Wrapper::IndexType;

  std::vector< std::vector<IndexType> > indices;
  std::vector< std::vector<double> > distances;

  Wrapper nnWrapper( container );
  nnWrapper.neighbourSearch( k, indices, distances );

  for( decltype(n) i = 0; i < n; i++ )
  {
    double density  = aleph::math::accumulate_kahan( distances[i].begin(), distances[i].end(), 0.0 );
    density         = -density;
    density        /= static_cast<double>( n );

    densities.push_back( density );
  }

  return densities;
}

} // namespace aleph

#endif
