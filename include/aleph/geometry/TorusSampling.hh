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
  Diaconis et al., samples exactly $n$ points from a torus with an inner
  radius of \f$R\f$ and an outer radius of \f$r\f$.

  @param R Inner radius
  @param r Outer radius
  @param n Number of samples to draw

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

  std::uniform_real_distribution<T> xDistribution  ( T(0), T(2.0 * M_PI) ); // distribution for x-values
  std::uniform_real_distribution<T> yDistribution  ( T(0), T(1.0 / M_PI) ); // distribution for y-values
  std::uniform_real_distribution<T> psiDistribution( T(0), T(2 * M_PI) );   // distribution for angles

  while( angles.size() < n )
  {
    auto x = xDistribution( rngX );
    auto y = yDistribution( rngY );
    auto f = static_cast<T>( 1.0 + (r/R) * std::cos( x ) ) / ( 2.0 * M_PI );

    if( y < f )
      angles.push_back( std::make_pair( x, psiDistribution( rngPsi ) ) );
  }

  return angles;
}

/**
  Converts a vector of angles into a point cloud that contains samples
  from a torus.

  @param angles Sampled angles to create the torus
  @param R Major radius
  @param r Minor radius

  @returns Point cloud containing the sampled points from the torus
*/

template <class T> aleph::containers::PointCloud<T> makeTorus( const std::vector< std::pair<T, T> >& angles, T R, T r )
{
  using PointCloud = aleph::containers::PointCloud<T>;

  PointCloud pc( angles.size(), 3 );

  unsigned index = 0;
  for( auto&& pair : angles )
  {
    // Note that this terminology follows Diaconis et al. who use "\psi"
    // insteadf of "\phi".
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
