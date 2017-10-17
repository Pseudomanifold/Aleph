#ifndef ALEPH_PERSISTENCE_DIAGRAMS_ENTROPY_HH__
#define ALEPH_PERSISTENCE_DIAGRAMS_ENTROPY_HH__

#include <aleph/containers/PointCloud.hh>

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

} // namespace aleph

#endif
