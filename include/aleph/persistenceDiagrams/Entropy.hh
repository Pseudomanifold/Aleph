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

/**
  Auxiliary function for calculating logarithms for a basis of `2`
  while handling a value of zero gracefully. This follows the way
  the logarithm is usually defined in statistics.

  @param x Input value
  @returns Logarithm of \p x to the basis of `2`
*/

template <class T> T log2( T x )
{
  if( x == T() )
    return T();
  else
    return std::log2( x );
}

/**
  @class RegularGrid
  @brief Auxiliary class for representing a regular grid
*/

template <class T> class RegularGrid
{
public:
  RegularGrid( unsigned width, unsigned height,
               T x0, T x1,
               T y0, T y1 )
    : _width( width )
    , _height( height )
    , _x0( x0 ), _x1( x1 ), _xOffset( (_x1 - _x0) / ( _width - 1 ) )
    , _y0( y0 ), _y1( y1 ), _yOffset( (_y1 - _y0) / ( _height - 1 ) )
    , _cells( new unsigned[ _width * _height ] )
  {
    std::fill( this->begin(), this->end(), 0 );
  }

  ~RegularGrid()
  {
    delete[] _cells;
  }

  // I do not want to implement these functions because they are
  // currently not required by the code below.
  RegularGrid( const RegularGrid& other )            = delete;
  RegularGrid& operator=( const RegularGrid& other ) = delete;

  unsigned* begin() { return _cells; }
  unsigned* end()   { return _cells + _width * _height; }

  unsigned& operator()( T x, T y )
  {
    x = x - _x0;
    y = y - _y0;

    unsigned i = unsigned( x / _xOffset );
    unsigned j = unsigned( y / _yOffset );

    return this->operator()(i,j);
  }

  unsigned& operator()( unsigned i, unsigned j )
  {
    return _cells[j * _width + i];
  }

  unsigned size() const noexcept
  {
    return _width * _height;
  }

private:
  unsigned _width;
  unsigned _height;

  T _x0, _x1, _xOffset;
  T _y0, _y1, _yOffset;

  unsigned* _cells;
};

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

/**
  Calculates a spatial entropy measure based on gridding data (or
  *quadrat counting*). Here, the idea is to measure the intensity
  of every grid cell, i.e. the number of points it contains. This
  quantity is then used as a probability. This probability can be
  used to define a notion of entropy subsequently.

  @param n Number of subdivisions of the grid in both directions
  @returns Grid-based spatial entropy
*/

template <class T> T gridEntropy( const aleph::PersistenceDiagram<T>& diagram, unsigned n )
{
  // Transform the data first in order to align the grid better with the
  // structure of the persistence points.
  auto pc = detail::makePointCloud( diagram, true );

  std::vector<T> X;
  std::vector<T> Y;

  X.reserve( pc.size() );
  Y.reserve( pc.size() );

  for( std::size_t i = 0; i < pc.size(); i++ )
  {
    auto&& p = pc[i];
    auto   x = p.front();
    auto   y = p.back();

    X.push_back(  0.5 * std::sqrt(2) * x + 0.5 * std::sqrt(2) * y );
    Y.push_back( -0.5 * std::sqrt(2) * x + 0.5 * std::sqrt(2) * y );
  }

  auto minmax_x = std::minmax_element( X.begin(), X.end() );
  auto minmax_y = std::minmax_element( Y.begin(), Y.end() );

  if( X.empty() || Y.empty() )
    return T();

  detail::RegularGrid<T> grid( n, n,
                               *minmax_x.first, *minmax_x.second,
                               *minmax_y.first, *minmax_y.second );

  for( std::size_t i = 0; i < X.size(); i++ )
    grid( X[i], Y[i] ) += 1;

  std::vector<T> entropies( grid.size() );

  std::transform( grid.begin(), grid.end(), entropies.begin(),
                  [&pc] ( unsigned n )
                  {
                    if( n != 0 )
                    {
                      T p = n / static_cast<T>( pc.size() );
                      T e = p * detail::log2( p );

                      return e;
                    }
                    else
                      return T();
                  } );

  return -aleph::math::accumulate_kahan_sorted( entropies.begin(),
                                                entropies.end(),
                                                T() );
}

} // namespace aleph

#endif
