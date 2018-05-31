/*
  This is a tool shipped by 'Aleph - A Library for Exploring Persistent
  Homology'.

  Its purpose is to create a *cover tree* for an input point cloud
  of arbitrary dimensionality. The cover tree creation process can
  be manipulated to some extent.
*/

#include <aleph/geometry/CoverTree.hh>
#include <aleph/geometry/Point.hh>

#include <aleph/geometry/distances/Euclidean.hh>
#include <aleph/geometry/distances/Wrapper.hh>

#include <aleph/utilities/String.hh>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <istream>
#include <stdexcept>
#include <string>
#include <vector>

#include <getopt.h>

using DataType  = double;
using Point     = aleph::geometry::Point<DataType>;
using Distance  = aleph::geometry::distances::Euclidean<DataType>;
using Wrapper   = aleph::geometry::distances::Wrapper<Distance, Point>;
using CoverTree = aleph::geometry::CoverTree<Point, Wrapper>;

std::vector<Point> load( std::istream& in )
{
  using namespace aleph::utilities;

  std::vector<Point> points;
  std::string line;

  while( std::getline( in, line ) )
  {
    line        = trim( line );
    auto tokens = splitByWhitespace( line );

    // Skip empty lines and comment lines
    if( line.empty() || line.front() == '#' || tokens.empty() )
      continue;

    std::vector<DataType> values;

    std::transform( tokens.begin(), tokens.end(),
      std::back_inserter( values ),
      [] ( const std::string& token )
      {
        return static_cast<DataType>( std::stod( token ) );
      }
    );

    points.emplace_back( Point( values.begin(), values.end() ) );
  }

  if( !points.empty() )
  {
    auto d = points.front().dimension();
    for( auto&& p : points )
      if( p.dimension() != d )
        throw std::runtime_error( "Dimension of points must not vary" );
  }

  return points;
}

std::vector<Point> load( const std::string& filename )
{
  std::ifstream in( filename );
  if( !in )
    throw std::runtime_error( "Unable to open input filename" );

  return load( in );
}

int main( int argc, char** argv )
{
  // Linkage criterion to use for constructing a hierarchical graph from
  // the cover tree. Currently, this is *not* used.
  std::string linkage = "single";

  {
    static option commandLineOptions[] =
    {
      { "linkage", required_argument, nullptr, 'l' },
      { nullptr  ,                 0, nullptr,  0  }
    };

    int option = 0;
    while( ( option = getopt_long( argc, argv, "l:", commandLineOptions, nullptr ) ) != -1 )
    {
      switch( option )
      {
      case 'l':
        linkage = optarg;
        break;
      default:
        throw std::runtime_error( "Unknown command-line option" );
      }
    }
  }

  // The remaining command-line arguments are considered to be file
  // names. The name '-' deserves special handling.
  std::string filename = argv[ optind++ ];

  std::vector<Point> points;
  if( filename == "-" )
    points = load( std::cin );
  else
    points = load( filename );

  // Nothing to do here
  if( points.empty() )
    return 0;

  std::cerr << "* Read "
            << points.size()
            << " points of dimension "
            << points.front().dimension()
            << "\n";

  CoverTree ct;

  ct.insert( points.begin(), points.end() );

  ct.print( std::cerr );

  for( auto&& p : points )
    ct.checkDistance( p );

  //for( auto&& p : points )
  //{
  //  ct.insert( p );

  //  if( !ct.isHarmonic( p ) )
  //    std::cerr << "Attempting to restore harmony\n";
  //}

  //ct.print( std::cerr );
}
