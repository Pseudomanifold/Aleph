#ifndef ALEPH_TOPOLOGY_IO_VTK_HH__
#define ALEPH_TOPOLOGY_IO_VTK_HH__

#include <cstddef>

#include <fstream>
#include <regex>
#include <stdexcept>
#include <sstream>
#include <string>
#include <vector>

#include "utilities/String.hh"

namespace aleph
{

namespace topology
{

namespace io
{

/**
  @class VTKStructuredGridReader
  @brief Simple reader class for VTK structured grids

  This class is a simple parser for VTK files in 'legacy format'. It is
  capable of parsing a structured grid and converting it to a simplicial
  complex. Data and weights of the simplicial complex will be taken from
  the VTK file.
*/

class VTKStructuredGridReader
{
public:
  template <class SimplicialComplex> void operator()( const std::string& filename, SimplicialComplex& K )
  {
    std::ifstream in( filename );
    if( !in )
      throw std::runtime_error( "Unable to read input file" );

    this->operator()( in, K );
  }

  template <class SimplicialComplex> void operator()( std::ifstream& in, SimplicialComplex& /* K */ )
  {
    using namespace aleph::utilities;

    using Simplex    = typename SimplicialComplex::ValueType;
    using DataType   = typename Simplex::DataType;
    using VertexType = typename Simplex::VertexType;

    std::string line;

    // Parse header first ----------------------------------------------

    std::size_t nx, ny, nz,
                n,
                s;

    bool parsedHeader = this->parseHeader( in, nx, ny, nz, n, s );
    if( !parsedHeader )
      return;

    // TODO: Check data type size against 's' and report if there are
    // issues such as insufficient storage space

    // Parse body ------------------------------------------------------

    std::vector<DataType> values;
    values.reserve( n );

    while( std::getline( in, line ) )
    {
      line         = trim( line );
      auto strings = split( line );

      for( auto&& value : strings )
        values.push_back( convert<DataType>( value ) );
    }

    VertexType v;
    (void) v;

    // Create topology -------------------------------------------------
    //
    // The first dimension (x) is increasing fastest. We first need some
    // mapping functions that assign indices based on x,y,z offsets.
    //
    // nx = 3
    // ny = 3
    // nz = 3
    //
    // [0,0,0] [1,0,0] [2,0,0] |  0, 1, 2
    // [0,1,0] [1,1,0] [2,1,0] |  3, 4, 5
    // [0,2,0] [1,2,0] [2,2,0] |  6, 7, 8
    // [0,0,1] [1,0,1] [2,0,1] |  9,10,11
    // [0,1,1] [1,1,1] [2,1,1] | 12,13,14
    // [0,2,1] [1,2,1] [2,2,1] | 15,16,17
    // [0,0,2] [1,0,2] [2,0,2] | 18,19,20
    // [0,1,2] [1,1,2] [2,1,2] | 21,22,23
    // [0,2,2] [1,2,2] [2,2,2] | 24,25,26
    //
    // nx = 1
    // ny = 2
    // nz = 3
    //
    // [0,0,0] [0,1,0] | 0,1
    // [0,0,1] [0,1,1] | 2,3
    // [0,0,2] [0,1,2] | 4,5
    //
    // [0,0,0] | 0
    // [0,1,0] | 1
    // [0,0,1] | 2
    // [0,1,1] | 3
    // [0,0,2] | 4
    // [0,1,2] | 5

#if 0
    auto indexToCoordinates = [nx,ny] ( std::size_t i, std::size_t& x, std::size_t& y, std::size_t& z )
    {
      x = i % nx;
      y = static_cast<std::size_t>( i / (nx)    ) % ny;
      z = static_cast<std::size_t>( i / (nx*ny) );
    };

    auto coordinatesToIndex = [nx,ny] ( std::size_t x, std::size_t y, std::size_t z )
    {
      return z * nx*ny + x % nx + y * nx;
    };
#endif

  }

private:

  /**
    Converts x,y,z coordinates to the corresponding index in the array
    of values.
  */

  static std::size_t coordinatesToIndex( const std::size_t nx, const std::size_t ny,
                                         const std::size_t x , const std::size_t  y, const std::size_t z ) noexcept
  {
   return z * nx*ny + x % nx + y * nx;
  }

  /**
    Enumerates all valid neighbours of a vertex and returns their
    indices. A valid neighbour is a neighbour whose index doesn't
    exceed the bounds of the grid. At present, the diagonal isn't
    used.
  */

  static std::vector<std::size_t> neighbours( const std::size_t nx, const std::size_t ny, const std::size_t nz,
                                              std::size_t x, std::size_t y, std::size_t z )
  {
    std::vector<std::size_t> neighbours;
    neighbours.reserve( 6 );

    // left
    if( x > 0 )
      neighbours.push_back( coordinatesToIndex(nx, ny, x-1, y, z) );
    // right
    else if( x+1 < nx )
      neighbours.push_back( coordinatesToIndex(nx, ny, x+1, y, z) );

    // bottom
    if( y > 0 )
      neighbours.push_back( coordinatesToIndex(nx, ny, x, y-1, z) );
    // top
    else if( y+1 < ny )
      neighbours.push_back( coordinatesToIndex(nx, ny, x, y+1, z) );

    // back
    if( z > 0 )
      neighbours.push_back( coordinatesToIndex(nx, ny, x, y, z-1) );
    // front
    else if( z+1 < nz )
      neighbours.push_back( coordinatesToIndex(nx, ny, x, y, z+1) );

    return neighbours;
  }

  /**
    Attempts parsing the header of a structured VTK file. If successful,
    returns true and sets all output variables:

    - x: x dimension
    - y: y dimension
    - z: z dimension
    - n: number of data points

    This function expects the header to conform to the VTK legacy file
    format.
  */

  bool parseHeader( std::ifstream& in,
                    std::size_t& x, std::size_t& y, std::size_t& z,
                    std::size_t& n,
                    std::size_t& s )
  {
    using namespace aleph::utilities;

    std::string identifier; // e.g. '# vtk DataFile Version 3.0'
    std::string header;     // e.g. 'VTK Output'
    std::string format;     // e.g. 'ASCII'
    std::string structure;  // e.g. 'DATASET STRUCTURED_GRID'
    std::string dimensions; // e.g. 'DIMENSIONS 100 10 1'
    std::string attributes; // e.g. 'POINT_DATA 1000 float'

    std::getline( in, identifier );
    std::getline( in, header );
    std::getline( in, format );

    if( !in )
      return false;

    identifier = trim( identifier );
    format     = trim( format );

    // This identifier is a little bit more lenient than the original
    // documentation requires: it will also accept if some fields are
    // joined by multiple spaces.
    std::regex reIdentifier( "#[[:space:]]+vtk[[:space:]]+DataFile[[:space:]]+Version[[:space:]]+([[:digit:]]+)\\.([[:digit:]]+)" );

    if( !std::regex_match( identifier, reIdentifier ) )
      return false;

    if( format != "ASCII" )
      throw std::runtime_error( "Binary file parsing is not yet supported" );

    std::getline( in, structure );
    std::getline( in, dimensions );
    std::getline( in, attributes );

    if( !in )
      return false;

    structure  = trim( structure );
    dimensions = trim( dimensions );
    attributes = trim( attributes );

    std::regex reStructure( "DATASET[[:space:]]+STRUCTURED_GRID" );
    if( !std::regex_match( structure, reStructure ) )
      return false;

    std::regex reDimensions( "DIMENSIONS[[:space:]]+([[:digit:]]+)[[:space:]]+([[:digit:]]+)[[:space:]]+([[:digit:]]+)" );
    std::smatch matches;

    if( !std::regex_match( dimensions, matches, reDimensions ) )
      return false;

    auto sx = matches[1];
    auto sy = matches[2];
    auto sz = matches[3];

    x = std::size_t( std::stoull( sx ) );
    y = std::size_t( std::stoull( sy ) );
    z = std::size_t( std::stoull( sz ) );

    std::regex reAttributes( "POINTS[[:space:]]+([[:digit:]]+)[[:space:]]+([[:alpha:]]+)" );
    if( !std::regex_match( attributes, matches, reAttributes ) )
      return false;

    auto sn = matches[1];
    auto st = matches[2];

    n = std::size_t( std::stoull( sn ) );

    if( st == "double" )
      s = sizeof(double);
    else if( st == "float" )
      s = sizeof(float);
    else if( st == "long" )
      s = sizeof(long);
    else if( st == "unsigned_long" )
      s = sizeof(unsigned long);
    else if( st == "int" )
      s = sizeof(int);
    else if( st == "unsigned_int" )
      s = sizeof(unsigned int);
    else if( st == "short" )
      s = sizeof(short);
    else if( st == "unsigned_short" )
      s = sizeof(unsigned short);
    else if( st == "char" )
      s = sizeof(char);
    else if( st == "unsigned char" )
      s = sizeof(unsigned char);
    else if( st == "bit" )
      s = sizeof(bool);

    return true;
  }
};

} // namespace io

} // namespace topology

} // namespace aleph

#endif
