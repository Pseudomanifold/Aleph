#ifndef ALEPH_GEOMETRY_SPHERE_SAMPLING_HH__
#define ALEPH_GEOMETRY_SPHERE_SAMPLING_HH__

#include <aleph/containers/PointCloud.hh>

#include <cmath>

#include <random>
#include <vector>

namespace aleph
{

namespace geometry
{

/**
  Samples \f$n\f$ points from a sphere such that the expected number of
  points per area of the sphere is uniform. Only the angular values of
  the sampled points (\f$\theta\f$, \f$\phi\f$) will be returned.

  @param n Number of samples to draw

  @returns Vector of angle values, i.e. \f$\theta\f$ and \f$\phi\f$,
           which are sufficient to describe the sphere. Please use
           aleph::geometry::makeSphere() to create  a point cloud
           from the resulting angles.
*/

template <class T>
std::vector< std::pair<T, T> > sphereSampling( unsigned n )
{
  std::random_device rd;
  std::mt19937 rng( rd() );

  std::vector< std::pair<T, T> > angles;
  angles.reserve( n );

  std::uniform_real_distribution<T> distribution( std::nextafter( T(0), std::numeric_limits<T>::max() ), T(1) );

  for( unsigned i = 0; i < n; i++ )
  {
    auto u = distribution( rng );
    auto v = distribution( rng );

    T theta = T( 2*M_PI*u );
    T phi   = std::acos( 1- 2*v );

    angles.push_back( std::make_pair( theta, phi ) );
  }

  return angles;
}

/**
  Converts a vector of angles into a point cloud that contains samples
  from a sphere of a given radius.

  @param angles The sampled angles of the sphere; these are the only
                data required to actually build the sphere, the other
                parameters merely control scaling

  @param r  Radius of the sphere
  @param x0 Centre x-position of the sphere
  @param y0 Centre y-position of the sphere
  @param z0 Centre z-position of the sphere

  @returns Point cloud containing the sphere samples
*/

template <class T> aleph::containers::PointCloud<T> makeSphere( const std::vector< std::pair<T, T> >& angles,
                                                                T r,
                                                                T x0 = T(),
                                                                T y0 = T(),
                                                                T z0 = T() )
{
  using PointCloud = aleph::containers::PointCloud<T>;

  PointCloud pc( angles.size(), 3 );

  unsigned index = 0;
  for( auto&& pair : angles )
  {
    auto theta = pair.first;
    auto phi   = pair.second;

    auto x = x0 + r * std::cos( theta ) * std::sin( phi );
    auto y = y0 + r * std::sin( theta ) * std::sin( phi );
    auto z = z0 + r * std::cos( phi );

    pc.set( index++, {x,y,z} );
  }

  return pc;
}

} // namespace geometry

} // namespace aleph

#endif
