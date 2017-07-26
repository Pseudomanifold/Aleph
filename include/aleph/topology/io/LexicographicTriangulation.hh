#ifndef ALEPH_TOPOLOGY_IO_LEXICOGRAPHIC_TRIANGULATION_HH__
#define ALEPH_TOPOLOGY_IO_LEXICOGRAPHIC_TRIANGULATION_HH__

#include <aleph/utilities/String.hh>

#include <algorithm>
#include <fstream>
#include <istream>
#include <regex>
#include <string>
#include <stdexcept>

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

  template <class SimplicialComplex, class OutputIterator> void operator()( const std::string& filename, OutputIterator result )
  {
    std::ifstream in( filename );
    if( !in )
      throw std::runtime_error( "Unable to read input file" );

    this->operator()( in, result );
  }

  template <class SimplicialComplex, class OutputIterator> void operator()( std::istream& in, OutputIterator result )
  {
    using namespace aleph::utilities;

    std::string line;
    std::string block; // contains _all_ data belonging to the current block

    Mode mode = Mode::ParsingBlocks;

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
        {
          // FIXME
          throw "Mode::ParsingBlocks";
        }

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
        // TODO: parse block
      }
    }

    (void) result;
  }

private:

  enum class Mode
  {
    ParsingBlocks,
    ParsingList
  };

  static bool isBlockFinished( const std::string& block ) noexcept
  {
    auto numOpeningBrackets = std::count( block.begin(), block.end(), '[' );
    auto numClosingBrackets = std::count( block.begin(), block.end(), ']' );

    return numOpeningBrackets == numClosingBrackets;
  }
};

} // namespace io

} // namespace topology

} // namespace aleph

#endif
