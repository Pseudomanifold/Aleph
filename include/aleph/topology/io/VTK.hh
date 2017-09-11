#ifndef ALEPH_TOPOLOGY_IO_VTK_HH__
#define ALEPH_TOPOLOGY_IO_VTK_HH__

#include <cstddef>

#include <algorithm>
#include <fstream>
#include <regex>
#include <stdexcept>
#include <sstream>
#include <string>
#include <vector>

#include <aleph/utilities/String.hh>

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
  the VTK file. Various query functions permit reading the name and data
  type of scalars, for example.
*/

class VTKStructuredGridReader
{
public:

  /**
    Reads a simplicial complex from a file, using the default maximum
    functor for weight assignment. If you want to change the functor,
    please refer to the overloaded variant of this method.

    @param filename Input filename
    @param  K       Simplicial complex

    @see operator()( const std::string&, SimplicialComplex& )
  */

  template <class SimplicialComplex> void operator()( const std::string& filename, SimplicialComplex& K )
  {
    using Simplex    = typename SimplicialComplex::ValueType;
    using DataType   = typename Simplex::DataType;

    this->operator()( filename, K, [] ( DataType a, DataType b ) { return std::max(a,b); } );
  }

  /**
    Reads a simplicial complex from a file while supporting arbitrary
    functors for weight assignment. The functor needs to support this
    interface:

    \code{.cpp}
    using SimplexType = typename SimplicialComplex::ValueType;
    using DataType    = typename Simplex::DataType;

    DataType Functor::operator()( DataType a, DataType b )
    {
      // Do something with a and b. Typically, this would be calculating
      // either the minimum or the maximum...
      return std::max(a, b);
    }
    \endcode

    Please refer to the documentation of SimplicialComplexReader::operator()( const std::string& SimplicialComplex&, Functor )
    for more details.
  */

  template <class SimplicialComplex, class Functor> void operator()( const std::string& filename, SimplicialComplex& K, Functor f )
  {
    std::ifstream in( filename );
    if( !in )
      throw std::runtime_error( "Unable to read input file" );

    this->operator()( in, K, f );
  }

  /** @overload operator()( const std::string&, SimplicialComplex& ) */
  template <class SimplicialComplex> void operator()( std::ifstream& in, SimplicialComplex& K )
  {
    using Simplex    = typename SimplicialComplex::ValueType;
    using DataType   = typename Simplex::DataType;

    this->operator()( in, K, [] ( DataType a, DataType b ) { return std::max(a,b); } );
  }

  /** @see operator()( const std::string&, SimplicialComplex&, SimplicialComplex&, Functor ) */
  template <class SimplicialComplex, class Functor> void operator()( std::ifstream& in, SimplicialComplex& K, Functor f )
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

    // This stores the data type size and makes it possible for a client
    // to look it up and compare it to the requested size of the complex
    // in order to see whether it is capable of storing the data.
    _dataTypeSize = s;

    // Parse body ------------------------------------------------------
    //
    // The body contains coordinates for each of the points, which are
    // dutifully ignored for now, and attributes. For now, point-based
    // attributes are supported.

    std::regex rePointData( "POINT_DATA[[:space:]]+([[:digit:]]+)" );
    std::regex reScalars( "SCALARS[[:space:]]+([[:alnum:]]+)[[:space:]]+([[:alnum:]]+)[[:space:]]*([[:digit:]]*)" );
    std::regex reLookupTable( "LOOKUP_TABLE[[:space:]]+([[:alnum:]]+)" );

    std::vector<DataType> coordinates;
    coordinates.reserve( n*3 );

    while( std::getline( in, line ) )
    {
      line         = trim( line );
      auto strings = split( line );

      for( auto&& coordinate : strings )
        coordinates.push_back( convert<DataType>( coordinate ) );

      if( coordinates.size() == n*3 )
        break;
    }

    std::vector<DataType> values;
    values.reserve( n );

    {
      std::smatch matches;

      while( std::getline( in, line ) )
      {
        if( std::regex_match( line, matches, rePointData ) )
        {
          auto sn = matches[1];
          if( static_cast<std::size_t>( std::stoull( sn ) ) != n )
            throw std::runtime_error( "Format error: number of point data attributes does not match number of points" );
        }
        else if( std::regex_match( line, matches, reScalars ) )
        {
          auto name       = matches[1];
          auto type       = matches[2];
          auto components = matches[3];

          _scalarsName = name;
          _scalarsType = type;

          if( !std::string( components ).empty() && std::stoul( components ) != 1 )
            throw std::runtime_error( "Format error: cannot handle scalars with more than one component" );
        }
        else if( std::regex_match( line, matches, reLookupTable ) )
        {
          auto name = matches[1];
          if( name != "default" )
            throw std::runtime_error( "Handling non-default lookup tables is not yet implemented" );
        }
        else
        {
          line         = trim( line );
          auto strings = split( line );

          for( auto&& value : strings )
            values.push_back( convert<DataType>( value ) );
        }
      }
    }

    // Create topology -------------------------------------------------
    //
    // Notice that this class only adds 0-simplices and 1-simplices to
    // the simplicial complex for now. While it is possible to include
    // triangles (i.e. 2-simplices), their creation order is not clear
    // and may subtly influence calculations.

    std::vector<Simplex> simplices;

    // Create 0-simplices ----------------------------------------------

    {
      VertexType v = VertexType();

      for( auto&& value : values )
        simplices.push_back( Simplex(v++, value) );
    }

    // Create 1-simplices ----------------------------------------------

    for( std::size_t z = 0; z < nz; z++ )
    {
      for( std::size_t y = 0; y < ny; y++ )
      {
        for( std::size_t x = 0; x < nx; x++ )
        {
          auto i = coordinatesToIndex(nx,ny,x,y,z);
          auto N = neighbours(nx,ny,nz,x,y,z);

          for( auto&& j : N )
          {
            if( j > i )
              continue;

            auto wi = simplices.at(i).data();
            auto wj = simplices.at(j).data();

            // Use the functor specified by the client in order to
            // assign a weight for the new simplex.
            auto w  = f(wi, wj);

            simplices.push_back( Simplex( {i,j}, w ) );
          }
        }
      }
    }

    K = SimplicialComplex( simplices.begin(), simplices.end() );
  }

  /** @returns Last read data type size */
  std::size_t dataTypeSize() const noexcept { return _dataTypeSize; }

  /** @returns Last read scalars name */
  std::string scalarsName() const noexcept  { return _scalarsName; }

  /** @returns Last read scalars data type */
  std::string scalarsType() const noexcept  { return _scalarsType; }

private:

  /**
    Converts an index in the array of values to the corresponding set of
    coordinates.
  */

  static void indexToCoordinates ( const std::size_t nx, const std::size_t ny,
                                   const std::size_t i,
                                   std::size_t& x,      // out: x
                                   std::size_t& y,      // out: y
                                   std::size_t& z )     // out: z
  {
    x = i % nx;
    y = static_cast<std::size_t>( i / (nx)    ) % ny;
    z = static_cast<std::size_t>( i / (nx*ny) );
  };

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
    if( x+1 < nx )
      neighbours.push_back( coordinatesToIndex(nx, ny, x+1, y, z) );

    // bottom
    if( y > 0 )
      neighbours.push_back( coordinatesToIndex(nx, ny, x, y-1, z) );
    // top
    if( y+1 < ny )
      neighbours.push_back( coordinatesToIndex(nx, ny, x, y+1, z) );

    // back
    if( z > 0 )
      neighbours.push_back( coordinatesToIndex(nx, ny, x, y, z-1) );
    // front
    if( z+1 < nz )
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
    std::string points;     // e.g. 'POINTS 1000 float'

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
    std::getline( in, points );

    if( !in )
      return false;

    structure  = trim( structure );
    dimensions = trim( dimensions );
    points     = trim( points );

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

    std::regex rePoints( "POINTS[[:space:]]+([[:digit:]]+)[[:space:]]+([[:alpha:]]+)" );
    if( !std::regex_match( points, matches, rePoints ) )
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

  /** Stores last read data type size */
  std::size_t _dataTypeSize = 0;

  /** Stores last read name of scalars (if present) */
  std::string _scalarsName;

  /** Stores last read data type of scalars (if present) */
  std::string _scalarsType;
};

} // namespace io

} // namespace topology

} // namespace aleph

#endif
