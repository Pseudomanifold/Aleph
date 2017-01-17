#ifndef ALEPH_CONTAINERS_DATA_DESCRIPTORS_HH__
#define ALEPH_CONTAINERS_DATA_DESCRIPTORS_HH__

#include <cmath>

#include <numeric>
#include <vector>

#include "distances/Euclidean.hh"
#include "distances/Traits.hh"

namespace aleph
{

template <class Distance, class Container> std::vector<double> eccentricities( const Container& container,
                                                                               unsigned order = 1 )
{
  auto n = container.size();
  auto d = container.dimension();

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
        distance = std::pow( distance, decltype(distance)( order ) ) / decltype(distance)(n);

        distances.push_back( distance );
      }
    }

    eccentricity = std::accumulate( distances.begin(), distances.end(), 0.0 );
    eccentricity = std::pow( eccentricity, 1.0 / order );

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
    auto p       = container[i];
    auto density = 0.0;

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

} // namespace aleph

#endif
