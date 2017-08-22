#ifndef ALEPH_GEOMETRY_TORUS_SAMPLING_HH__
#define ALEPH_GEOMETRY_TORUS_SAMPLING_HH__

#include <aleph/containers/PointCloud.hh>

#include <random>
#include <vector>

#include <cmath>

namespace aleph
{

namespace geometry
{

/**
  Using the rejection sampling method from "Sampling from a manifold" by
  Diaconis et al., samples at most $n$ points from a torus with an inner
  radius of \f$R\f$ and an outer radius of \f$r\f$.

  @param R Inner radius
  @param r Outer radius
  @param n Maximum number of samples to draw

  @returns Vector of angle values, i.e. \f$\theta\f$ and \f$\psi\f$,
           which are sufficient to describe a torus. Please use
           aleph::geometry::makeTorus() to create a point cloud
           from the resulting angles.
*/

template <class T>
std::vector< std::pair<T, T> > torusRejectionSampling( T R,
                                                       T r,
                                                       unsigned n )
{
  std::random_device rd;

  std::mt19937 rngPsi( rd() );
  std::mt19937 rngX( rd() );
  std::mt19937 rngY( rd() );

  // I do not store the values of theta directly, but instead report directly
  // those angles that are deemed to be "correct".
  std::vector< std::pair<T, T> > angles;
  angles.reserve( n );

  std::vector<T> xValues;
  std::vector<T> yValues;
  xValues.reserve( n );
  yValues.reserve( n );

  std::uniform_real_distribution<T> xDistribution( T(0), T(2.0 * M_PI) );
  std::uniform_real_distribution<T> yDistribution( T(0), T(1.0 / M_PI) );

  for( unsigned i = 0; i < n; i++ )
  {
    xValues.push_back( xDistribution( rngX ) );
    yValues.push_back( yDistribution( rngY ) );
  }

  std::vector<double> functionValues;
  functionValues.reserve( n );

  for( double x : xValues )
    functionValues.push_back( static_cast<T>( 1.0 + (r/R) * std::cos( x ) ) / ( 2.0 * M_PI ) );

  std::uniform_real_distribution<T> psiDistribution( T(0), T(2 * M_PI) );

  for( unsigned i = 0; i < n; i++ )
    if( yValues[i] < functionValues[i] )
      angles.push_back( std::make_pair( xValues[i], psiDistribution( rngPsi ) ) );

  return angles;

}

/**
  Converts a vector of angles into a point cloud that contains samples
  from a torus.

  @param angles Sampled angles to create the torus
  @param R Inner radius
  @param r Outer radius

  @returns Point cloud containing the sampled points from the torus
*/

template <class T> aleph::containers::PointCloud<T> makeTorus( const std::vector< std::pair<T, T> >& angles, T R, T r )
{
  using PointCloud = aleph::containers::PointCloud<T>;

  PointCloud pc( angles.size(), 3 );

  unsigned index = 0;
  for( auto&& pair : angles )
  {
    auto theta = pair.first;
    auto psi   = pair.second;

    auto x = ( R + r * std::cos( theta ) ) * std::cos( psi );
    auto y = ( R + r * std::cos( theta ) ) * std::sin( psi );
    auto z = (     r * std::sin( theta ) );

    pc.set( index++, {x,y,z} );
  }

  return pc;
}

} // namespace geometry

} // namespace aleph

#endif
