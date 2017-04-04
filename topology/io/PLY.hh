#ifndef ALEPH_TOPOLOGY_IO_PLY_HH__
#define ALEPH_TOPOLOGY_IO_PLY_HH__

#include "filtrations/Data.hh"
#include "utilities/String.hh"

#include "topology/Simplex.hh"
#include "topology/SimplicialComplex.hh"

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

namespace topology
{

namespace io
{

namespace detail
{

/*
  Maps PLY data types to their corresponding sizes in bytes. This is
  required for parsing binary files later on.
*/

std::map<std::string, unsigned short> TypeSizeMap = {
  { "double" , 8 },
  { "float"  , 4 },
  { "int"    , 4 },
  { "int32"  , 4 },
  { "uint"   , 4 },
  { "uint32" , 4 },
  { "short"  , 2 },
  { "ushort" , 2 },
  { "char"   , 1 },
  { "uchar"  , 1 },
  { "uint8"  , 1 }
};

} // namespace detail

/**
  @class PLYReader
  @brief Parses PLY files

  This is a simple reader class for files in PLY format. It supports
  reading PLY files with an arbitrary number of vertex properties. A
  user may specify which property to use in order to assign the data
  stored for each simplex.
*/

class PLYReader
{
public:

  // Contains all descriptors for a single PLY property. The number of
  // elements is somewhat superfluous when parsing ASCII files.
  struct PropertyDescriptor
  {
    unsigned index; // Offset of attribute
    unsigned bytes; // Number of bytes
  };

  template <class SimplicialComplex> void operator()( const std::string& filename, SimplicialComplex& K )
  {
    std::ifstream in( filename );
    if( !in )
      throw std::runtime_error( "Unable to read input file" );

    this->operator()( in, K );
  }

  template <class SimplicialComplex> void operator()( std::ifstream& in, SimplicialComplex& K )
  {
    // Header ------------------------------------------------------------
    //
    // The header needs to consist of the word "ply", followed by a "format"
    // description.

    std::size_t numVertices = 0;
    std::size_t numFaces    = 0;

    // Current line in file. This is required because I prefer reading the
    // file line by line via `std::getline`.
    std::string line;

    bool headerParsed       = false;
    bool parseBinary        = false;
    bool littleEndian       = false;

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

      if( format == "ascii 1.0" )
        parseBinary = false;
      else if( format == "binary_little_endian 1.0" )
      {
        parseBinary  = true;
        littleEndian = true;
      }
      else if( format == "binary_big_endian 1.0" )
      {
        parseBinary  = true;
        littleEndian = false;
      }
      else
        throw std::runtime_error( "Format error: Expecting \"ascii 1.0\" or \"binary_little_endian 1.0\" or \"binary_big_endian 1.0\" " );
    }

    // Maps properties of a PLY file to an index. The index specifies at
    // which position in a single vertex specification line the selected
    // property appears.
    //
    // The parser expects certain properties, viz. "x", "y", and "z" to be
    // present in all files. Else, an error is raised.
    std::map<std::string, PropertyDescriptor> properties;
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

        // Skip property lists; we will handle them implicitly when
        // converting all faces.
        if( dataType == "list" )
          continue;

        if( !converter )
          throw std::runtime_error( "Property conversion error: Expecting data type and name of property" );

        name = utilities::trim( name );

        PropertyDescriptor descriptor;
        descriptor.index = propertyIndex;
        descriptor.bytes = detail::TypeSizeMap.at( dataType );

        properties[name] = descriptor;

        ++propertyIndex;
      }

      if( line == "end_header" )
      {
        headerParsed = true;
        break;
      }
    }
    while( !headerParsed && in );

    assert( numVertices > 0 );
    assert( numFaces    > 0 );

    using Simplex = typename SimplicialComplex::ValueType;

    // Container for storing all simplices that are created while reading
    // the mesh data structure.
    std::vector<Simplex> simplices;

    if( parseBinary )
    {
    }
    else
      simplices = this->parseASCII<Simplex>( in, numVertices, numFaces, properties );

    in.close();

    K = SimplicialComplex( simplices.begin(), simplices.end() );
    K.recalculateWeights();
    K.sort( filtrations::Data<Simplex>() );
  }

  /* Sets the property to read for every simplex */
  void setDataProperty( const std::string& property )
  {
    _property = property;
  }

private:

  template <class Simplex> std::vector<Simplex> parseASCII( std::ifstream& in,
                                                            std::size_t numVertices, std::size_t numFaces,
                                                            const std::map<std::string,PropertyDescriptor>& properties )
  {
    using DataType   = typename Simplex::DataType;
    using VertexType = typename Simplex::VertexType;

    std::vector<Simplex> simplices;
    std::string line;

    // Read vertices -----------------------------------------------------

    for( std::size_t vertexIndex = 0; vertexIndex < numVertices; vertexIndex++ )
    {
      std::vector<double> vertexCoordinates( 3 );

      std::getline( in, line );

      line        = utilities::trim( line );
      auto tokens = utilities::split( line );

      auto ix     = properties.at( "x" ).index;
      auto iy     = properties.at( "y" ).index;
      auto iz     = properties.at( "z" ).index;

      auto x      = std::stod( tokens.at( ix ) );
      auto y      = std::stod( tokens.at( iy ) );
      auto z      = std::stod( tokens.at( iz ) );

      _coordinates.push_back( {x,y,z} );

      // No property for reading weights specified, or the specified
      // property could not be found; just use the default weight of
      // the simplex class.
      if( _property.empty() || properties.find( _property ) == properties.end() )
        simplices.push_back( { VertexType( vertexIndex ) } );
      else
      {
        auto iw    = properties.at( _property ).index;
        DataType w = aleph::utilities::convert<DataType>( tokens.at(iw) );

        simplices.push_back( Simplex( VertexType( vertexIndex ), w ) );
      }
    }

    // Read faces --------------------------------------------------------

    // Keep track of all edges that are encountered. This ensures that the
    // simplicial complex is valid upon construction and does not have any
    // missing simplices.
    std::set< std::pair<VertexType, VertexType> > edges;

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

    return simplices;
  }

  /** Data property to assign to new simplices */
  std::string _property = "z";

  /** Coordinates stored by the current run of the reader */
  std::vector< std::vector<double> > _coordinates;
};

} // namespace io

} // namespace topology

} // namespace aleph

#endif
