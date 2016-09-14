#ifndef ALEPH_PLY_HH__
#define ALEPH_PLY_HH__

#include "Simplex.hh"
#include "SimplicialComplex.hh"

#include "filtrations/Data.hh"
#include "utilities/String.hh"

#include <cassert>

#include <fstream>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace aleph
{

namespace io
{

template <
  class DataType,
  class VertexType
> SimplicialComplex< Simplex<DataType, VertexType> > loadPLY( const std::string& filename, const std::string& property = std::string() )
{
  std::ifstream in( filename );

  if( !in )
    throw std::runtime_error( "Unable to open input filename" );

  using Simplex           = Simplex<DataType, VertexType>;
  using SimplicialComplex = SimplicialComplex<Simplex>;

  // This stores the coordinates for further processing. They are not required
  // as the simplicial complex only uses the connectivity information. However,
  // the coordinates are of course helpful when rendering the complex.
  std::vector< std::vector<double> > coordinates;

  // Container for storing all simplices that are created while reading
  // the mesh data structure.
  std::vector<Simplex> simplices;

  // Keep track of all edges that are encountered. This ensures that the
  // simplicial complex is valid upon construction and does not have any
  // missing simplices.
  std::set< std::pair<VertexType, VertexType> > edges;

  // Current line in file. This is required because I prefer reading the
  // file line by line via `std::getline`.
  std::string line;

  // Header ------------------------------------------------------------
  //
  // The header needs to consist of the word "ply", followed by a "format"
  // description.

  std::size_t numVertices = 0;
  std::size_t numFaces    = 0;

  bool headerParsed       = false;

  std::getline( in, line );
  line = utilities::trim( line );

  if( line != "ply" )
    throw std::runtime_error( "Format error: Expecting \"ply\"" );

  std::getline( in, line );
  line = utilities::trim( line );

  if( line.substr( 0, 6 ) != "format" )
    throw std::runtime_error( "Format error: Expecting \"format\"" );
  else
  {
    std::string format = line.substr( 6 );
    format = utilities::trim( format );

    if( format != "ascii 1.0" )
      throw std::runtime_error( "Format error: Expecting \"ascii 1.0\"" );
  }

  // Maps properties of a PLY file to an index. The index specifies at
  // which position in a single vertex specification line the selected
  // property appears.
  //
  // The parser expects certain properties, viz. "x", "y", and "z" to be
  // present in all files. Else, an error is raised.
  std::map<std::string, unsigned> propertyMap;

  unsigned propertyIndex = 0;

  // Parse the rest of the header, taking care to skip any comment lines.
  do
  {
    std::getline( in, line );
    line = utilities::trim( line );

    if( !in )
      break;

    if( line.substr( 0, 7 ) == "comment" )
      continue;
    else if( line.substr( 0, 7) == "element" )
    {
      std::string element = line.substr( 7 );
      element = utilities::trim( element );

      std::istringstream converter( element );

      std::string name;
      std::size_t numElements = 0;

      converter >> name
                >> numElements;

      if( !converter )
        throw std::runtime_error( "Element conversion error: Expecting number of elements" );

      name = utilities::trim( name );

      if( name == "vertex" )
        numVertices = numElements;
      else if( name == "face" )
        numFaces = numElements;
    }
    else if( line.substr( 0, 8) == "property" )
    {
      std::string property = line.substr( 8 );
      property = utilities::trim( property );

      std::istringstream converter( property );

      std::string dataType;
      std::string name;

      converter >> dataType
                >> name;

      if( !converter )
        throw std::runtime_error( "Property conversion error: Expecting data type and name of property" );

      name                = utilities::trim( name );
      propertyMap[ name ] = propertyIndex;

      ++propertyIndex;
    }

    if( line == "end_header" )
    {
      headerParsed = true;
      break;
    }
  }
  while( !headerParsed );

  assert( numVertices > 0 );
  assert( numFaces    > 0 );

  // Read vertices -----------------------------------------------------

  for( std::size_t vertexIndex = 0; vertexIndex < numVertices; vertexIndex++ )
  {
    std::vector<double> vertexCoordinates( 3 );

    std::getline( in, line );

    line        = utilities::trim( line );
    auto tokens = utilities::split( line );

    auto ix     = propertyMap.at( "x" );
    auto iy     = propertyMap.at( "y" );
    auto iz     = propertyMap.at( "z" );

    auto x      = std::stod( tokens.at( ix ) );
    auto y      = std::stod( tokens.at( iy ) );
    auto z      = std::stod( tokens.at( iz ) );

    coordinates.push_back( {x,y,z} );

    // No property for reading weights specified; just add a simplex
    // with the default weight.
    if( property.empty() )
      simplices.push_back( { VertexType( vertexIndex ) } );
    else
    {
      // FIXME: Check for existence first
      // FIXME: Conversion is stupid
      auto iw    = propertyMap.at( property );
      DataType w = static_cast<DataType>( std::stod( tokens.at( iw ) ) );

      simplices.push_back( Simplex( VertexType( vertexIndex ), w ) );
    }
  }

  // Read faces --------------------------------------------------------

  for( std::size_t faceIndex = 0; faceIndex < numFaces; faceIndex++ )
  {
    std::getline( in, line );
    std::istringstream converter( line );

    unsigned numEntries = 0;

    converter >> numEntries;
    if( !converter )
      throw std::runtime_error( "Face conversion error: Expecting number of entries" );

    // I can make a simplex out of a triangle, but every other shape would
    // get complicated.
    if( numEntries == 3 )
    {
      VertexType i1 = 0;
      VertexType i2 = 0;
      VertexType i3 = 0;

      converter >> i1
                >> i2
                >> i3;

      if( !converter )
        throw std::runtime_error( "Unable to parse vertex indices" );

      Simplex triangle( {i1,i2,i3} );

      // Create edges ----------------------------------------------

      for( auto itEdge = triangle.begin_boundary();
           itEdge != triangle.end_boundary();
           ++itEdge )
      {
        // As the boundary iterator works as a filtered iterator only,
        // I need this copy.
        //
        // Else, different calls to `begin()` and `end()` will result in
        // two different copies of the simplex. The copies, in turn,
        // will then not be usable as a source for vertices.
        Simplex edge = *itEdge;

        auto u = *( edge.begin() );
        auto v = *( edge.begin() + 1 );

        if( u < v )
          std::swap( u, v );

        auto pair = edges.insert( std::make_pair( u,v ) );

        if( pair.second )
          simplices.push_back( Simplex( edge.begin(), edge.end() ) );
      }

      simplices.push_back( triangle );
    }
    else
      throw std::runtime_error( "Format error: Expecting triangular faces only" );
  }

  in.close();

  SimplicialComplex K( simplices.begin(), simplices.end() );
  K.recalculateWeights();
  K.sort( filtrations::Data<Simplex>() );

  return K;
}

}

}

#endif
