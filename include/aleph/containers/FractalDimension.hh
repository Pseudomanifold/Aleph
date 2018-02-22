#ifndef ALEPH_CONTAINERS_FRACTAL_DIMENSION_HH__
#define ALEPH_CONTAINERS_FRACTAL_DIMENSION_HH__

#include <aleph/geometry/distances/Traits.hh>

#include <aleph/math/Statistics.hh>

#include <algorithm>
#include <map>
#include <stdexcept>
#include <vector>

#include <cmath>

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
  auto n = container.size();
  auto d = container.dimension();

  std::map<double, unsigned> distances;

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

      distances[ traits.from( distance ) ] += 1;
    }
  }

  CorrelationDimensionSequence cds;
  cds.x.reserve( distances.size() );
  cds.y.reserve( distances.size() );

  // Determine the correlation dimension integral for all potential
  // values. This only requires counting how many *pairs* have been
  // seen by the algorithm.

  unsigned seen = 0;
  for( auto&& pair : distances )
  {
    auto&& distance = pair.first;
    auto&& count    = pair.second;

    seen += count;

    if( distance > 0 )
    {
      cds.x.push_back( distance );
      cds.y.push_back( seen / ( static_cast<double>( n * (n-1) ) * 0.5 ) );
    }
  }

  return cds;
}

/**
  Estimates the correlation dimension from a correlation dimension
  sequence, which involves calculating a log-log plot of the data,
  and determining the best coefficient for a linear fit.
*/

double correlationDimension( const CorrelationDimensionSequence& cds )
{
  if( cds.x.size() != cds.y.size() )
    throw std::runtime_error( "Inconsistent correlation dimension sequence" );

  std::vector<double> X;
  std::vector<double> Y;

  X.reserve( cds.x.size() );
  Y.reserve( cds.y.size() );

  auto log = [] ( double x )
  {
    return std::log( x );
  };

  std::transform( cds.x.begin(), cds.x.end(), std::back_inserter( X ), log );
  std::transform( cds.y.begin(), cds.y.end(), std::back_inserter( Y ), log );

  // This is a simple linear regression model. We are only interested in
  // the slope of the regression line, so this is sufficient.

  auto cov
    = aleph::math::sampleCovariance( X.begin(), X.end(),
                                     Y.begin(), Y.end() );

  auto var
    = aleph::math::sampleVariance( X.begin(), X.end() );

  return cov / var;
}

} // namespace containers

} // namespace aleph

#endif
