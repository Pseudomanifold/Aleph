/*
  This is an example file shipped by 'Aleph - A Library for Exploring
  Persistent Homology'.

  This example demonstrates how to create random persistence diagrams
  that may be used to compare topological algorithms with each other,
  such as persistence indicator functions and Wasserstein distances.

  Demonstrated classes:

    - aleph::distances::Euclidean
    - aleph::geometry::BruteForce
    - aleph::PersistenceDiagram

  Demonstrated functions:

    - aleph::calculatePersistenceDiagrams
    - aleph::geometry::buildVietorisRipsComplex
    - aleph::geometry::makeSphere
    - aleph::geometry::makeTorus
    - aleph::geometry::sphereSampling
    - aleph::geometry::torusRejectionSampling

  Original author: Bastian Rieck
*/

#include <aleph/distances/Euclidean.hh>

#include <aleph/geometry/BruteForce.hh>
#include <aleph/geometry/SphereSampling.hh>
#include <aleph/geometry/TorusSampling.hh>
#include <aleph/geometry/VietorisRipsComplex.hh>

#include <aleph/persistenceDiagrams/PersistenceDiagram.hh>

#include <aleph/persistentHomology/Calculation.hh>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <limits>
#include <random>
#include <sstream>

#include <cmath>

#include <getopt.h>

// We first have to specify the data type of the persistence diagram,
// i.e. the type that is used by its individual points.
using DataType           = double;
using Distance           = aleph::distances::Euclidean<DataType>;
using PointCloud         = aleph::containers::PointCloud<DataType>;
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

/**
  Auxiliary function for creating a random persistence diagram based on
  random samples from a torus. The torus is described by two radii, one
  for the outer part, the other for the inner part.

  Due to the sampling technique used, the specified number is merely an
  upper bound for the number of points that are to be sampled. Moreover
  the function will automatically handle Vietoris--Rips expansion.
*/

PersistenceDiagram createRandomTorusPersistenceDiagram( DataType R, DataType r, unsigned n )
{
  auto pointCloud = aleph::geometry::makeTorus(
    aleph::geometry::torusRejectionSampling( R, r, n ),
    R, r
  );

  aleph::geometry::BruteForce<PointCloud, Distance> bruteForceWrapper( pointCloud );

  auto K
    = aleph::geometry::buildVietorisRipsComplex( bruteForceWrapper, r, 2 );

  auto diagrams
    = aleph::calculatePersistenceDiagrams( K );

  diagrams.at(1).removeDiagonal();

  // We are only interested in the one-dimensional persistent homology
  // of the samples.
  return diagrams.at(1);
}

/**
  Auxiliary function for creating a random persistence diagram based on
  random samples from a sphere with a given radius.

  Note that this function automatically handles Vietoris--Rips
  expansion.
*/

PersistenceDiagram createRandomSpherePersistenceDiagram( DataType r, unsigned n )
{
  auto pointCloud = aleph::geometry::makeSphere(
    aleph::geometry::sphereSampling<DataType>( n ),
    r
  );

  aleph::geometry::BruteForce<PointCloud, Distance> bruteForceWrapper( pointCloud );

  auto K
    = aleph::geometry::buildVietorisRipsComplex( bruteForceWrapper, r, 2 );

  auto diagrams
    = aleph::calculatePersistenceDiagrams( K );

  diagrams.at(1).removeDiagonal();

  // We are only interested in the one-dimensional persistent homology
  // of the samples.
  return diagrams.at(1);
}

int main( int argc, char** argv )
{
  static option commandLineOptions[] = {
      { "m"     , required_argument, nullptr, 'm' },
      { "n"     , required_argument, nullptr, 'n' },
      { "R"     , required_argument, nullptr, 'R' },
      { "r"     , required_argument, nullptr, 'r' },
      { "sphere", no_argument      , nullptr, 's' },
      { "torus" , no_argument      , nullptr, 't' },
      { nullptr , 0                , nullptr,  0  }
  };

  unsigned m           = 50;
  unsigned n           = 50;
  DataType R           = DataType(0.25);
  DataType r           = DataType(0.50);

  bool sampleFromSphere = false;
  bool sampleFromTorus  = false;

  int option = 0;
  while( ( option = getopt_long( argc, argv, "m:n:R:r:t", commandLineOptions, nullptr ) ) != -1 )
  {
    switch( option )
    {
    case 'm':
      m = static_cast<unsigned>( std::stoul(optarg) );
      break;
    case 'n':
      n = static_cast<unsigned>( std::stoul(optarg) );
      break;
    case 'R':
      R = static_cast<DataType>( std::stod(optarg) );
      break;
    case 'r':
      r = static_cast<DataType>( std::stod(optarg) );
      break;
    case 's':
      sampleFromSphere = true;
      sampleFromTorus  = false;
      break;
    case 't':
      sampleFromSphere = false;
      sampleFromTorus  = true;
      break;
    default:
      break;
    }
  }

  std::cerr << "* Sampling " << n << " persistence diagrams\n";
  if( sampleFromSphere )
    std::cerr << "* Sampling " << m << " points from a sphere with r=" << r << "\n";
  else if( sampleFromTorus )
    std::cerr << "* Sampling at most " << m << " points from a torus with R=" << R << " and r=" << r << "\n";
  else
    std::cerr << "* Generating " << m << " random points per diagram\n";

  for( unsigned i = 0; i < n; i++ )
  {
    PersistenceDiagram pd;

    if( sampleFromSphere )
      pd = createRandomSpherePersistenceDiagram(r, m);
    else if( sampleFromTorus )
      pd = createRandomTorusPersistenceDiagram(R, r, m);
    else
      pd = createRandomPersistenceDiagram(m);

    std::stringstream stream;
    stream << "/tmp/";
    stream << std::setw( static_cast<int>( std::floor( std::log10(n)+1 ) ) ) << std::setfill( '0' ) << i;
    stream << ".txt";

    std::ofstream out( stream.str() );
    out << pd << "\n";
  }
}
