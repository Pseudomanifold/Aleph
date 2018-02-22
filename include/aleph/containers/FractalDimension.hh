#ifndef ALEPH_CONTAINERS_FRACTAL_DIMENSION_HH__
#define ALEPH_CONTAINERS_FRACTAL_DIMENSION_HH__

#include <aleph/geometry/distances/Traits.hh>

#include <vector>

namespace aleph
{

namespace containers
{

/*
 Wrapper class for representing a sequence of correlation dimension
 integral values. I want the interface to be as clear as possible.
*/

struct CorrelationDimensionSequence
{
  std::vector<double> x;
  std::vector<double> y;
};

/**
  Calculates samples of the correlation dimension integral for a given
  point cloud. This does *not* yet result in a dimension estimate, but
  only produces a set of points.
*/

template <
  class Distance,
  class Container
> CorrelationDimensionSequence correlationDimensionIntegral(
  const Container& container,
  Distance dist = Distance() )
{
  aleph::geometry::distances::Traits<Distance> traits;

  auto n = container.size();
  auto d = container.dimension();

  std::vector<double> distances;
  distances.reserve( n * (n - 1) / 2 );

  aleph::geometry::distances::Traits<Distance> traits;

  for( decltype(n) i = 0; i < n; i++ )
  {
    auto&& p = container[i];

    for( decltype(n) j = i+1; j < n; j++ )
    {
      auto&& q      = container[j];
      auto distance = dist( p.begin(),
                            q.begin(),
                            d );

      distances.push_back( trait.from( distance ) );
    }
  }
}

} // namespace containers

} // namespace aleph

#endif
