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

#include <aleph/persistenceDiagrams/Entropy.hh>
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
using Distance           = aleph::geometry::distances::Euclidean<DataType>;

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
    , _x0( x0 ), _x1( x1 ), _xOffset( (_x1 - _x0) / ( _width - 1 ) )
    , _y0( y0 ), _y1( y1 ), _yOffset( (_y1 - _y0) / ( _height - 1 ) )
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

DataType gridEntropy( const PointCloud& pc, unsigned n )
{
  std::vector<DataType> X;
  std::vector<DataType> Y;

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
    return DataType();

  RegularGrid grid( n, n,
                    *minmax_x.first, *minmax_x.second,
                    *minmax_y.first, *minmax_y.second );

  for( std::size_t i = 0; i < X.size(); i++ )
    grid( X[i], Y[i] ) += 1;

  std::vector<DataType> entropies( grid.size() );

  std::transform( grid.begin(), grid.end(), entropies.begin(),
                  [&pc] ( unsigned n )
                  {
                    if( n != 0 )
                    {
                      DataType p = n / static_cast<DataType>( pc.size() );
                      DataType e = p * log( p );

                      return e;
                    }
                    else
                      return DataType();
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
  {
    auto e_nn = aleph::nearestNeighbourAreaEntropy( input.persistenceDiagram );
    auto e_rg = gridEntropy( input.pointCloud, 20 );

    std::cout << e_nn << "\t" << e_rg << "\n";
  }
}
