#ifndef ALEPH_TOPOLOGY_IO_LEXICOGRAPHIC_TRIANGULATION_HH__
#define ALEPH_TOPOLOGY_IO_LEXICOGRAPHIC_TRIANGULATION_HH__

#include <aleph/utilities/String.hh>

#include <algorithm>
#include <fstream>
#include <istream>
#include <iterator>
#include <regex>
#include <string>
#include <stdexcept>
#include <vector>

#include <cctype>

namespace aleph
{

namespace topology
{

namespace io
{

/**
  @class LexicographicTriangulationReader
  @brief Class for reading a list of lexicographic triangulations

  This class supports loading triangulations in the lexicographic format
  developed by Frank H. Lutz. The format contains an identifier for each
  manifold, followed by a list of simplices. Start and end of every item
  is denoted by "[" and "]", respectively.

  An example of a valid triangulation:

  \code
  manifold_2_4_1=[[1,2,3],[1,2,4],[1,3,4],[2,3,4]]
  \endcode

  This parser supports whitespace characters at any point of the list or
  the simplex. Moreover, comments indicated by '#' at the beginning of a
  line are supported as well.

  @see http://page.math.tu-berlin.de/~lutz/stellar/mixed.html
*/

class LexicographicTriangulationReader
{
public:

  template <class SimplicialComplex> void operator()( const std::string& filename, std::vector<SimplicialComplex>& result )
  {
    std::ifstream in( filename );
    if( !in )
      throw std::runtime_error( "Unable to read input file" );

    this->operator()( in, result );
  }

  template <class SimplicialComplex> void operator()( std::istream& in, std::vector<SimplicialComplex>& result )
  {
    using namespace aleph::utilities;

    std::string line;
    std::string block; // contains _all_ data belonging to the current block

    Mode mode = Mode::ParsingBlocks;

    result.clear();
    result.shrink_to_fit();

    while( std::getline( in, line ) )
    {
      line     = trim( line );
      auto pos = line.find( '=' );

      // Skip empty lines and comment lines
      if( line.empty() || line.front() == '#' )
        continue;

      if( pos != std::string::npos )
      {
        block.clear();

        if( mode != Mode::ParsingBlocks )
          throw std::runtime_error( "Format error; unexpected open block detected" );

        mode  = Mode::ParsingList;
        block = block + line.substr(pos+1);
        block = trim( block );
      }
      else if( mode == Mode::ParsingList )
      {
        block = block + line;
        block = trim( block );
      }

      if( isBlockFinished( block ) )
      {
        result.emplace_back( parseBlock<SimplicialComplex>( block ) );
        mode      = Mode::ParsingBlocks;

        block.clear();
      }
    }

    result.shrink_to_fit();
  }

private:

  enum class Mode
  {
    ParsingBlocks,  // Parsing blocks of triangulations
    ParsingList     // Parsing lists of simplices
  };

  static bool isBlockFinished( const std::string& block ) noexcept
  {
    auto numOpeningBrackets = std::count( block.begin(), block.end(), '[' );
    auto numClosingBrackets = std::count( block.begin(), block.end(), ']' );

    return numOpeningBrackets == numClosingBrackets;
  }

  /**
    Parses a block consisting of a list of simplices. Since the block
    cannot be easily tokenized at the beginning, this function merely
    uses a simple lookahead method.

    Upon success, a non-empty simplicial complex is returned.
  */

  template <class SimplicialComplex> static SimplicialComplex parseBlock( const std::string& block ) noexcept
  {
    // Sanity check
    if( !isBlockFinished( block ) )
      return {};

    using Simplex    = typename SimplicialComplex::ValueType;
    using VertexType = typename Simplex::VertexType;

    SimplicialComplex K;

    std::vector<VertexType> vertices; // list of vertices belonging to the current simplex
    bool parsingSimplex = false;      // flag indicating whether the parser is between two
                                      // simplices or parsing an individual simplex

    using DifferenceType = std::string::difference_type;
    using SizeType       = std::string::size_type;

    // Remove the opening bracket from the block. It makes the number of
    // brackets unbalanced but it also simplifies parsing because we may
    // rely on the fact that every opening bracket starts a new simplex.
    for( auto it = block.begin() + DifferenceType( block.find_first_of( '[' ) + 1 ); it != block.end(); )
    {
      // Either the beginning of a list of simplices, or the very
      // beginning of the whole block
      if( *it == '[' )
      {
        // Format error: the parser is already parsing a simplex
        if( parsingSimplex )
          return {};

        // The parser is now in simplex parsing mode
        else
          parsingSimplex = true;
      }

      // Either the end of a list of simplices, or the very end of the
      // whole block
      else if( *it == ']' )
      {
        // Add parsed simplex to the simplicial complex and clear the
        // set of vertices afterwards. Reset the parse mode as well.
        if( parsingSimplex )
        {
          K.push_back( Simplex( vertices.begin(), vertices.end() ) );
          vertices.clear();

          ++it;
          parsingSimplex = false;
        }
        // Finished parsing the simplicial complex
        else
          break;
      }

      // Just skip all valid characters and try to directly proceed with
      // the parsing process.
      if( *it == ',' || std::isspace( *it ) )
        ++it;

      // Look ahead until the next "]" appears. Tokenize the list of
      // vertices in between.
      if( parsingSimplex )
      {
        auto offset        = static_cast<SizeType>( std::distance( block.begin(), it ) );
        auto positionBegin = block.find_first_of( '[', offset );
        auto positionEnd   = block.find_first_of( ']', offset );

        // Format error: the list of simplices requires an opening and
        // a closing bracket
        if( positionBegin == std::string::npos || positionEnd == std::string::npos  )
          return {};

        using namespace aleph::utilities;

        auto list   = trim( block.substr( positionBegin + 1, positionEnd - positionBegin - 1 ) );
        auto tokens = split( list, std::string(",") );

        for( auto&& token : tokens )
          vertices.emplace_back( convert<VertexType>( token ) );

        std::advance( it, long( positionEnd - offset ) );
      }
    }

    return K;
  }
};

} // namespace io

} // namespace topology

} // namespace aleph

#endif
