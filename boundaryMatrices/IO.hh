#ifndef ALEPH_IO_HH__
#define ALEPH_IO_HH__

#include "BoundaryMatrix.hh"

#include <algorithm>
#include <fstream>
#include <iterator>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace aleph
{

template <class Representation> BoundaryMatrix<Representation> load( const std::string& filename )
{
  using Index = typename Representation::Index;

  std::ifstream in( filename );
  if( !in )
    return {};

  std::string line;
  Index numColumns = Index(0);

  while( std::getline( in, line ) )
  {
    // Remove all spaces
    line.erase( std::remove_if( line.begin(), line.end(), ::isspace ),
                line.end() );

    if( line.empty() == false && line.front() != '#' )
      ++numColumns;
  }

  in.clear();
  in.seekg( 0 );

  BoundaryMatrix<Representation> M;
  M.setNumColumns( numColumns );

  Index curColumn = Index(0);

  while( std::getline( in, line ) )
  {
    // Ignore empty lines and comment lines
    if( line.empty() || line.front() == '#' )
      continue;

    std::istringstream iss( line );

    std::vector<std::string> tokens;

    std::copy( std::istream_iterator<std::string>( iss ),
               std::istream_iterator<std::string>(),
               std::back_inserter( tokens ) );

    // Ignore this line if it does not contain any tokens
    if( tokens.empty() )
      continue;

    std::vector<Index> indices;
    indices.reserve( tokens.size() );

    std::transform( tokens.begin(), tokens.end(), std::back_inserter( indices ),
                    [] ( const std::string& token )
                    {
                      return static_cast<Index>( std::stoul( token ) );
                    } );

    if( indices.empty() )
      throw std::runtime_error( "Amount of indices in boundary must not be empty" );

    // TODO: This ignores the dimension of the column
    // TODO: This assumes that the column indices are ordered
    M.setColumn( curColumn,
                 indices.begin() + 1, indices.end() );

    if( M.getDimension( curColumn ) != indices.front() )
      throw std::runtime_error( "Inconsistency between actual number of indices and specified number of indices in boundary" );

    ++curColumn;
  }

  return M;
}

}

#endif
