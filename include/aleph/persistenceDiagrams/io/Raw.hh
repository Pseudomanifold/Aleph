#ifndef ALEPH_PERSISTENCE_DIAGRAMS_IO_RAW_HH__
#define ALEPH_PERSISTENCE_DIAGRAMS_IO_RAW_HH__

#include <aleph/persistenceDiagrams/PersistenceDiagram.hh>
#include <aleph/utilities/String.hh>

#include <fstream>
#include <string>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace aleph
{

namespace io
{

/**
  Loads a persistence diagram from a file. The file format is kept
  simple: lines starting with '#' are ignored, empty lines will be
  ignored as well.

  The function looks for two numeric tokens, separated by tabs, or
  by spaces, and adds them to the persistence diagram in the order
  in which they appear.

  Any errors will result in exceptions.
*/

template <class T> PersistenceDiagram<T> load( const std::string& filename )
{
  std::ifstream in( filename );

  if( !in )
    throw std::runtime_error( "Unable to open filename for reading" );

  using namespace aleph::utilities;

  PersistenceDiagram<T> persistenceDiagram;

  std::string line;
  while( std::getline( in, line ) )
  {
    line = trim( line );

    if( line.empty() || line.front() == '#' )
      continue;

    std::stringstream converter( line );

    std::string a_token;
    std::string b_token;

    converter >> a_token
              >> b_token;

    if( converter.fail() )
      throw std::runtime_error( "Unable to parse token" );
    else
    {
      T a = convert<T>( a_token );
      T b = convert<T>( b_token );

      persistenceDiagram.add( a, b );
    }
  }

  return persistenceDiagram;
}

} // namespace io

} // namespace aleph

#endif
