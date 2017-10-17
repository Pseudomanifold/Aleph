#ifndef ALEPH_PERSISTENCE_DIAGRAMS_ENTROPY_HH__
#define ALEPH_PERSISTENCE_DIAGRAMS_ENTROPY_HH__

#include <aleph/containers/PointCloud.hh>

#include <aleph/geometry/BruteForce.hh>
#include <aleph/geometry/FLANN.hh>

#include <aleph/geometry/distances/Euclidean.hh>

#include <aleph/persistenceDiagrams/PersistenceDiagram.hh>

#include <aleph/math/KahanSummation.hh>

#include <algorithm>
#include <vector>

#include <cmath>

namespace aleph
{

namespace detail
{

/**
  Converts a persistence diagram to a point cloud of the same data type.
  The diagram is optionally transformed into another *coordinate system*
  as suggested by Edelsbrunner et al. in the paper *Current Open Problems
  in Discrete and Computational Geometry*.

  @param diagram   Input persistence diagram
  @param transform Flag indicating whether the input shall be transformed

  @see http://www.mathnet.ru/eng/mais259
*/

template <class T> aleph::containers::PointCloud<T> makePointCloud( const aleph::PersistenceDiagram<T>& diagram, bool transform = false )
{
  aleph::containers::PointCloud<T> pc( diagram.size(), 2 );
  std::size_t i = 0;

  for( auto&& point : diagram )
  {
    std::vector<T> p;

    if( transform )
      p = { point.x() + point.y(), point.y() - point.x() };
    else
      p = { point.x(), point.y() };

    pc.set(i++, p.begin(), p.end() );
  }

  return pc;
}

template <class T> T log2( T x )
{
  if( x == T() )
    return T();
  else
    return std::log2( x );
}

} // namespace detail

/**
  Calculates the persistent entropy of a given persistence diagram. This
  notion of entropy was developed by Chintakunta et al. in the paper *An
  entropy-based persistence barcode*, Pattern Recognition Volume 48, No.
  2, pp. 391--401.

  @param D Persistence diagram
  @returns Persistent entropy

  @see https://doi.org/10.1016/j.patcog.2014.06.023
*/

template <class T> T persistentEntropy( const aleph::PersistenceDiagram<T>& D )
{
  using Point = typename aleph::PersistenceDiagram<T>::Point;

  aleph::math::KahanSummation<T> result = T();

  std::vector<T> persistenceValues;
  persistenceValues.reserve( D.size() );

  std::transform( D.begin(), D.end(),
                  std::back_inserter( persistenceValues ),
                  [] ( const Point& p )
                  {
                    return p.persistence();
                  } );

  auto totalPersistence = aleph::math::accumulate_kahan_sorted( persistenceValues.begin(), persistenceValues.end() );

  std::vector<T> probabilities;
  probabilities.reserve( D.size() );

  std::transform( persistenceValues.begin(), persistenceValues.end(),
                  std::back_inserter( probabilities ),
                  [&totalPersistence] ( T persistence )
                  {
                    auto p = persistence / totalPersistence;
                    return p * std::log2( p );
                  } );
}

/**
  Calculates a spatial entropy measure based on the distance to the
  nearest neighbour of every point in the diagram. This is based on
  the idea that a *large* distance between neighbours is consistent
  with a well-ordered diagram, whereas smaller distances indicate a
  clustering or clumping process.

  @param diagram Input diagram
  @returns Nearest-neighbour entropy
*/

template <class T> T nearestNeighbourAreaEntropy( const aleph::PersistenceDiagram<T>& diagram )
{
  using PointCloud = aleph::containers::PointCloud<T>;
  using Distance   = aleph::geometry::distances::Euclidean<T>;

#ifdef ALEPH_WITH_FLANN
  using NearestNeighbours = aleph::geometry::FLANN<PointCloud, Distance>;
#else
  using NearestNeighbours = aleph::geometry:::BruteForce<PointCloud, Distance>;
#endif

  auto pc = detail::makePointCloud( diagram );
  NearestNeighbours nn( pc );

  using ElementType = typename NearestNeighbours::ElementType;
  using IndexType   = typename NearestNeighbours::IndexType;

  std::vector< std::vector<IndexType> > indices;
  std::vector< std::vector<ElementType> > distances;

  nn.neighbourSearch( 2,  // the 1st nearest neighbour is the point itself
                      indices,
                      distances );

  std::vector<ElementType> nearestNeigbourDistances;
  nearestNeigbourDistances.reserve( pc.size() );

  for( auto&& distance : distances )
    nearestNeigbourDistances.push_back( distance.at(1) );

  std::vector<ElementType> areas( pc.size() );

  std::transform( nearestNeigbourDistances.begin(), nearestNeigbourDistances.end(),
                  areas.begin(),
                  [] ( ElementType radius )
                  {
                    return static_cast<ElementType>( radius * radius * 2 * M_PI );
                  } );

  auto totalArea
    = aleph::math::accumulate_kahan_sorted( areas.begin(), areas.end(),
                                            ElementType() );

  std::vector<T> entropies( pc.size() );

  std::transform( areas.begin(), areas.end(),
                  entropies.begin(),
                  [&totalArea] ( ElementType area )
                  {
                    T p = static_cast<T>( area / totalArea );
                    T e = p * detail::log2( p );

                    return e;
                  } );

  return -aleph::math::accumulate_kahan_sorted( entropies.begin(),
                                                entropies.end(),
                                                T() );

}

} // namespace aleph

#endif
