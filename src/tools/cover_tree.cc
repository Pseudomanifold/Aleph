/*
  This is a tool shipped by 'Aleph - A Library for Exploring Persistent
  Homology'.

  Its purpose is to create a *cover tree* for an input point cloud
  of arbitrary dimensionality. The cover tree creation process can
  be manipulated to some extent.
*/

#include <aleph/geometry/CoverTree.hh>
#include <aleph/geometry/Point.hh>

#include <aleph/utilities/String.hh>

#include <algorithm>
#include <fstream>
#include <istream>
#include <stdexcept>
#include <string>
#include <vector>

#include <getopt.h>

using DataType = double;
using Point    = aleph::geometry::Point<DataType>;

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

  return points;
}

std::vector<Point> load( const std::string& filename )
{
  std::ifstream in( filename );
  if( !in )
    throw std::runtime_error( "Unable to open input filename" );

  return load( in );
}



int main( int, char** )
{
}
