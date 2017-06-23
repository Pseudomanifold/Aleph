#ifndef ALEPH_TOPOLOGY_IO_PLY_HH__
#define ALEPH_TOPOLOGY_IO_PLY_HH__

#include <aleph/utilities/String.hh>

#include <aleph/topology/Simplex.hh>
#include <aleph/topology/SimplicialComplex.hh>

#include <aleph/topology/filtrations/Data.hh>

#include <cassert>

#include <fstream>
#include <limits>
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

/*
  Maps PLY data types to their 'signedness'. This is required for
  parsing binary files later on.

  TODO: Use this to extend parse
*/

std::map<std::string, bool> TypeSignednessMap = {
  { "double" , true },
  { "float"  , true },
  { "int"    , true },
  { "int32"  , true },
  { "uint"   , false },
  { "uint32" , false },
  { "short"  , true },
  { "ushort" , false },
  { "char"   , true },
  { "uchar"  , false },
  { "uint8"  , false }
};

/* Describes an arbitrary value of a PLY file */
union PLYValue
{
  double         d;
  float          f;
  int            i;
  short          s;
  char           c;
  unsigned       u;
  unsigned short us;
  unsigned char  uc;
};

/*
  Reads a single value from a binary input stream, reversing the storage
  order if necessary.
*/

template <class T> void readValue( std::ifstream& stream,
                                   std::size_t bytes,
                                   bool littleEndian,
                                   T* target )
{
  if( !littleEndian || bytes == 1 )
    stream.read( reinterpret_cast<char*>( target ), static_cast<std::streamsize>( bytes ) );
  else
  {
    char* buffer         = new char[bytes];
    char* reversedBuffer = new char[bytes];

    stream.read( buffer, static_cast<std::streamsize>( bytes ) );

    for( std::size_t i = 0; i < bytes; i++ )
      reversedBuffer[i] = buffer[ bytes - 1 - i ];

    std::copy( reversedBuffer, reversedBuffer + bytes, target );

    // FIXME: Superfluous?
#if 0
    memcpy( reinterpret_cast<void*>( &target ),
            reinterpret_cast<void*>( reversedBuffer ),
            bytes );
#endif

    delete[] reversedBuffer;
    delete[] buffer;
  }
}

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
    std::string name;     // Property name (or list name)
    unsigned index;       // Offset of attribute for ASCII data
    unsigned bytesOffset; // Offset of attribute for binary data
    unsigned bytes;       // Number of bytes

    // Only used for lists: Here, both the length parameter and the
    // entry parameter of a list usually have different lengths.
    unsigned bytesListSize;
    unsigned bytesListEntry;
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

    std::size_t numVertices       = 0; // Number of edges
    std::size_t vertexSizeInBytes = 0; // Only relevant for binary files
    std::size_t numFaces          = 0; // Number of faces
    std::size_t faceSizeInBytes   = 0; // Only relevant for binary files

    // Current line in file. This is required because I prefer reading the
    // file line by line via `std::getline`.
    std::string line;

    bool headerParsed = false;
    bool parseBinary  = false;
    bool littleEndian = false;

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

    // All properties stored in the PLY file in the order in which they
    // were discovered. The parse requires the existence of some of the
    // properties, e.g. "x", "y", and "z" in order to work correctly.
    std::vector<PropertyDescriptor> properties;

    unsigned propertyIndex  = 0; // Offset for properties in ASCII files
    unsigned propertyOffset = 0; // Offset for properties in binary files

    bool readingVertexProperties = false;
    bool readingFaceProperties   = false;

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
        {
          numVertices             = numElements;
          readingVertexProperties = true;
        }
        else if( name == "face" )
        {
          numFaces              = numElements;
          readingFaceProperties = true;
        }
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

        dataType = utilities::trim( dataType );
        name     = utilities::trim( name );

        PropertyDescriptor descriptor;
        descriptor.index = propertyIndex;

        // List of properties require a special handling. The syntax is
        // "property list SIZE_TYPE ENTRY_TYPE NAME", e.g. "property
        // list uint float vertex_height".
        if( dataType == "list" )
        {
          std::string sizeType = name;
          std::string entryType;
          std::string listName;

          converter >> entryType
                    >> listName;

          utilities::trim( entryType );
          utilities::trim( listName );

          descriptor.bytesListSize  = detail::TypeSizeMap.at( sizeType );
          descriptor.bytesListEntry = detail::TypeSizeMap.at( entryType );
          descriptor.name           = listName;
        }
        else
        {
          descriptor.bytes       = detail::TypeSizeMap.at( dataType );
          descriptor.bytesOffset = propertyOffset;
          descriptor.name        = name;
        }

        if( !converter )
          throw std::runtime_error( "Property conversion error: Expecting data type and name of property" );

        if( readingFaceProperties )
          faceSizeInBytes += descriptor.bytes;
        else if( readingVertexProperties )
          vertexSizeInBytes += descriptor.bytes;

        propertyOffset += descriptor.bytes;
        propertyIndex  += 1;

        properties.push_back( descriptor );
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
      simplices = this->parseBinary<Simplex>( in,
                                              numVertices, numFaces,
                                              littleEndian,
                                              properties );
    }
    else
    {
      simplices = this->parseASCII<Simplex>( in,
                                             numVertices, numFaces,
                                             properties );
    }

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

  template <class Simplex> std::vector<Simplex> parseBinary( std::ifstream& in,
                                                             std::size_t numVertices,
                                                             std::size_t numFaces,
                                                             bool littleEndian,
                                                             const std::vector<PropertyDescriptor>& properties )
  {
    using DataType   = typename Simplex::DataType;
    using VertexType = typename Simplex::VertexType;

    DataType weight = DataType();

    for( std::size_t vertexIndex = 0; vertexIndex < numVertices; vertexIndex++ )
    {
      std::vector<double> coordinates;

      for( auto&& descriptor : properties )
      {
        // Only faces may have lists for now...
        if( descriptor.bytesListSize + descriptor.bytesListEntry != 0 )
          continue;

        detail::PLYValue pv;

        switch( descriptor.bytes )
        {
          case 8:
            detail::readValue( in, descriptor.bytes, littleEndian, &pv.d );
            break;
          case 4:
            detail::readValue( in, descriptor.bytes, littleEndian, &pv.f );
            break;
          case 2:
            detail::readValue( in, descriptor.bytes, littleEndian, &pv.s );
            break;
          case 1:
            detail::readValue( in, descriptor.bytes, littleEndian, &pv.c );
            break;
        }

        // FIXME: Not sure whether this is legal...
        if( descriptor.name == "x" || descriptor.name == "y" || descriptor.name == "z" )
          coordinates.push_back( pv.d );
      }
    }

    for( std::size_t faceIndex = 0; faceIndex < numVertices; faceIndex++ )
    {
    }
  }


  template <class Simplex> std::vector<Simplex> parseASCII( std::ifstream& in,
                                                            std::size_t numVertices, std::size_t numFaces,
                                                            const std::vector<PropertyDescriptor>& properties )
  {
    using DataType   = typename Simplex::DataType;
    using VertexType = typename Simplex::VertexType;

    std::vector<Simplex> simplices;
    std::string line;

    auto getPropertyIndex = [&properties] ( const std::string& property )
    {
      auto it = std::find_if( properties.begin(), properties.end(),
                              [&property] ( const PropertyDescriptor& descriptor )
                              {
                                return descriptor.name == property;
                              } );

      if( it != properties.end() )
        return it->index;
      else
        return std::numeric_limits<unsigned>::max();
    };

    // Read vertices -----------------------------------------------------

    for( std::size_t vertexIndex = 0; vertexIndex < numVertices; vertexIndex++ )
    {
      std::vector<double> vertexCoordinates( 3 );

      std::getline( in, line );

      line        = utilities::trim( line );
      auto tokens = utilities::split( line );

      auto ix     = getPropertyIndex( "x" );
      auto iy     = getPropertyIndex( "y" );
      auto iz     = getPropertyIndex( "z" );
      auto iw     = getPropertyIndex( _property );

      auto x      = std::stod( tokens.at( ix ) );
      auto y      = std::stod( tokens.at( iy ) );
      auto z      = std::stod( tokens.at( iz ) );

      _coordinates.push_back( {x,y,z} );

      // No property for reading weights specified, or the specified
      // property could not be found; just use the default weight of
      // the simplex class.
      if( _property.empty() || iw == std::numeric_limits<unsigned>::max() )
        simplices.push_back( { VertexType( vertexIndex ) } );
      else
      {
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
