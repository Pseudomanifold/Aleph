#ifndef ALEPH_PERSISTENCE_DIAGRAMS_IO_JSON_HH__
#define ALEPH_PERSISTENCE_DIAGRAMS_IO_JSON_HH__

#include <aleph/config/RapidJSON.hh>

#include <aleph/persistenceDiagrams/PersistenceDiagram.hh>

#include <aleph/utilities/String.hh>

#include <fstream>
#include <istream>
#include <ostream>
#include <string>

#include <map>

#ifdef ALEPH_WITH_RAPID_JSON
  #include <rapidjson/document.h>
  #include <rapidjson/istreamwrapper.h>
#endif

namespace aleph
{

namespace io
{

/**
  Writes a persistence diagram to an output stream, using the JSON
  format. The diagram will be serialized such that its points will
  be stored in the field ``data`` as a two-dimensional array. Note
  that infinite values will be encoded as strings. Additional data
  about the diagram, e.g. its dimension, are stored in name--value
  pairs.

  An optional map can be supplied in order to store arbitrary data
  about each diagram.
*/

template <class Diagram> void writeJSON(
  std::ostream& o,
  const Diagram& D,
  const std::string& name = std::string(),
  const std::map<std::string, std::string>& kvs = std::map<std::string, std::string>() )
{
  std::string level = "  ";

  o << "{\n";

  o << level << "\"betti\": "     << D.betti()     << ",\n"
    << level << "\"dimension\": " << D.dimension() << ",\n";

  // Store additional key--value pairs belonging to the current diagram,
  // if they have been supplied by the client.
  if( !kvs.empty() )
  {
    for( auto&& pair : kvs )
      o << level << "\"" << pair.first << "\": " << "\"" << pair.second << "\",\n";
  }

  if( !name.empty() )
    o << level << "\"name\": " << "\"" << name << "\",\n";

  o << level << "\"size\": "      << D.size()      << ",\n"
    << level << "\"diagram\": "   << "[\n";

  for( auto it = D.begin(); it != D.end(); ++it )
  {
    if( it != D.begin() )
      o << ",\n";

    o << level << level << "["
                        << "\"" << it->x() << "\""
                        << ","
                        << "\"" << it->y() << "\""
                        << "]";
  }

  o << "\n"
    << level << "]\n"
    << "}";
}

/**
  Convenience function for writing a persistence diagram to a file
  in JSON format.
*/

template <class Diagram> void writeJSON( const std::string& filename, const Diagram& D, const std::string& name = std::string() )
{
  std::ofstream out( filename );
  if( !out )
    throw std::runtime_error( "Unable to open output file" );

  writeJSON( out, D, name );
}

/**
  Reads multiple persistence diagrams from an input stream in JSON
  format. The stream is checked for consistency; appropriate error
  messages will be raised if necessary.
*/

template <class T> std::vector< aleph::PersistenceDiagram<T> > readJSON( std::istream& in )
{
#ifdef ALEPH_WITH_RAPID_JSON

  rapidjson::Document document;

  {
    rapidjson::IStreamWrapper isw( in );
    document.ParseStream( isw );
  }

  using PersistenceDiagram = aleph::PersistenceDiagram<T>;

  std::vector<PersistenceDiagram> persistenceDiagrams;

  if( !document.HasMember( "diagrams" ) )
    throw std::runtime_error( "Unable to find array of persistence diagrams" );

  for( auto&& diagram : document["diagrams"].GetArray() )
  {
    auto dimension = diagram["dimension"].GetUint();
    auto betti     = diagram["betti"].GetUint();
    auto size      = diagram["size"].GetUint();

    PersistenceDiagram persistenceDiagram;
    persistenceDiagram.setDimension( dimension );

    for( auto&& point : diagram["diagram"].GetArray() )
    {
      auto x  = aleph::utilities::convert<T>( point[0].GetString() );
      auto y  = aleph::utilities::convert<T>( point[1].GetString() );

      persistenceDiagram.add( x,y );
    }

    if( persistenceDiagram.size() != size )
      throw std::runtime_error( "Stored number of points does not match number of points in persistence diagram" );

    if( persistenceDiagram.betti() != betti )
      throw std::runtime_error( "Stored Betti number does not match Betti number of persistence diagram" );

    persistenceDiagrams.emplace_back( persistenceDiagram );
  }

  return persistenceDiagrams;

#else

  // I'm only doing this because the compiler should not complain about
  // unused parameters.
  (void) in;

  return {};
#endif
}

/**
  Reads multiple persistence diagrams from a JSON input file. This
  is only a convenience unction provided to call the other reading
  routine defined above.
*/

template <class T> std::vector< aleph::PersistenceDiagram<T> > readJSON( const std::string& filename )
{
  std::ifstream in( filename );
  if( !in )
    throw std::runtime_error( "Unable to read input file" );

  return readJSON<T>( in );
}

} // namespace io

} // namespace aleph

#endif
