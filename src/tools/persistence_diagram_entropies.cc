/*
  This is a tool shipped by 'Aleph - A Library for Exploring Persistent
  Homology'.

  It permits the calculation of different entropy measures defined for
  persistence diagrams.

  Original author: Bastian Rieck
*/

#include <aleph/config/FLANN.hh>

#ifdef ALEPH_WITH_FLANN
  #include <aleph/geometry/FLANN.hh>
#endif

#include <aleph/containers/PointCloud.hh>

#include <aleph/geometry/BruteForce.hh>

#include <aleph/geometry/distances/Euclidean.hh>

#include <aleph/math/KahanSummation.hh>

#include <aleph/persistenceDiagrams/PersistenceDiagram.hh>

#include <aleph/persistenceDiagrams/io/Raw.hh>

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

#include <cmath>

using DataType           = double;
using PersistenceDiagram = aleph::PersistenceDiagram<DataType>;
using PointCloud         = aleph::containers::PointCloud<DataType>;
using Distance           = aleph::distances::Euclidean<DataType>;

#ifdef ALEPH_WITH_FLANN
  using NearestNeighbours = aleph::geometry::FLANN<PointCloud, Distance>;
#else
  using NearestNeighbours = aleph::geometry::BruteForce<PointCloud, Distance>;
#endif

class RegularGrid
{
public:
  RegularGrid( unsigned width, unsigned height,
               DataType x0, DataType x1,
               DataType y0, DataType y1 )
    : _width( width )
    , _height( height )
    , _x0( x0 ), _x1( x1 ), _xOffset( (_x1 - _x0) / _width )
    , _y0( y0 ), _y1( y1 ), _yOffset( (_y1 - _y0) / _height )
    , _cells( new unsigned[ _width * _height ] )
  {
    std::fill( this->begin(), this->end(), 0 );
  }

  unsigned* begin() { return _cells; }
  unsigned* end()   { return _cells + _width * _height; }

  unsigned& operator()( DataType x, DataType y )
  {
    x = x - _x0;
    y = y - _y0;

    unsigned i = unsigned( x / _xOffset + 0.5 * _xOffset );
    unsigned j = unsigned( y / _yOffset + 0.5 * _yOffset );

    if( i >= _width )
      i = _width - 1;

    if( j >= _height )
      j = _height - 1;

    return this->operator()(i,j);
  }

  unsigned& operator()( unsigned i, unsigned j )
  {
    return _cells[j * _width + i];
  }

private:
  unsigned _width;
  unsigned _height;

  DataType _x0, _x1, _xOffset;
  DataType _y0, _y1, _yOffset;

  unsigned* _cells;
};

struct Input
{
  std::string filename;
  PersistenceDiagram persistenceDiagram;
  PointCloud pointCloud;
};

PointCloud makePointCloud( const PersistenceDiagram& diagram )
{
  PointCloud pc( diagram.size(), 2 );
  std::size_t i = 0;

  for( auto&& point : diagram )
  {
    std::vector<DataType> p = { point.x(), point.y() };
    pc.set(i++, p.begin(), p.end() );
  }

  return pc;
}

template <class T> T log( T x )
{
  if( x == T() )
    return T();
  else
    return std::log( x );
}

DataType nearestNeighbourAreaEntropy( const PointCloud& pc )
{
  NearestNeighbours nn( pc );

  using ElementType = NearestNeighbours::ElementType;
  using IndexType   = NearestNeighbours::IndexType;

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

  std::vector<DataType> entropies( pc.size() );

  std::transform( areas.begin(), areas.end(),
                  entropies.begin(),
                  [&totalArea] ( ElementType area )
                  {
                    DataType p = static_cast<ElementType>( area / totalArea );
                    DataType e = p * log( p );

                    return e;
                  } );

  return -aleph::math::accumulate_kahan_sorted( entropies.begin(),
                                                entropies.end(),
                                                DataType() );
}

void usage()
{
}

int main( int argc, char** argv )
{
  if( argc <= 1 )
  {
    usage();
    return -1;
  }

  std::vector<Input> inputs;
  inputs.reserve( argc );

  for( int i = 1; i < argc; i++ )
  {
    std::string filename = argv[i];

    std::cerr << "* Loading '" << filename << "'...";

    auto&& diagram = aleph::io::load<DataType>( filename );

    Input input = {
      filename,
      diagram,
      makePointCloud( diagram )
    };

    inputs.push_back( input );

    std::cerr << "finished\n";
  }

  for( auto&& input : inputs )
    std::cout << nearestNeighbourAreaEntropy( input.pointCloud ) << "\n";
}
