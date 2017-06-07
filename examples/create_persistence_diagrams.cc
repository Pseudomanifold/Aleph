/*
  This is an example file shipped by 'Aleph - A Library for Exploring
  Persistent Homology'.

  This example demonstrates how to create random persistence diagrams
  that may be used to compare topological algorithms with each other,
  such as persistence indicator functions and Wasserstein distances.

  Demonstrated classes:

    - aleph::PersistenceDiagram

  Original author: Bastian Rieck
*/

#include "persistenceDiagrams/PersistenceDiagram.hh"

#include <iomanip>
#include <fstream>
#include <limits>
#include <random>
#include <sstream>

#include <cmath>

// We first have to specify the data type of the persistence diagram,
// i.e. the type that is used by its individual points.
using DataType           = double;
using PersistenceDiagram = aleph::PersistenceDiagram<DataType>;

/**
  Auxiliary function for creating a random persistence diagram. Points
  are drawn from a uniform distribution on [0,1]. The function ensures
  that all points are situated _above_ the diagonal.
*/

PersistenceDiagram createRandomPersistenceDiagram( unsigned n )
{
  std::random_device rd;
  std::default_random_engine rng( rd() );
  std::uniform_real_distribution<DataType> distribution( DataType(0), DataType( std::nextafter( DataType(1), std::numeric_limits<DataType>::max() ) ) );

  PersistenceDiagram D;

  for( unsigned i = 0; i < n; i++ )
  {
    auto x = distribution( rng );
    auto y = distribution( rng );

    if( x > y )
      std::swap( x,y );

    D.add( x,y );
  }

  return D;
}

int main( int, char** )
{
  // TODO: make configurable
  unsigned n = 50;
  unsigned m = 50;

  for( unsigned i = 0; i < n; i++ )
  {
    auto&& pd = createRandomPersistenceDiagram(m);

    std::stringstream stream;
    stream << "/tmp/";
    stream << std::setw( static_cast<int>( std::floor( std::log10(n)+1 ) ) ) << std::setfill( '0' ) << i;
    stream << ".txt";

    std::ofstream out( stream.str() );
    out << pd << "\n";
  }
}
