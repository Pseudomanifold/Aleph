#ifndef ALEPH_CONTAINERS_DATA_DESCRIPTORS_HH__
#define ALEPH_CONTAINERS_DATA_DESCRIPTORS_HH__

#include <cmath>

#include <numeric>
#include <vector>

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

}

#endif
