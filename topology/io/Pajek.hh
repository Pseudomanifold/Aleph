#ifndef ALEPH_TOPOLOGY_IO_PAJEK_HH__
#define ALEPH_TOPOLOGY_IO_PAJEK_HH__

#include <algorithm>
#include <fstream>
#include <regex>
#include <stdexcept>
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
  @class PajekReader
  @brief Parses files in Pajek format

  This is a simple reader for graphs in Pajek format. It supports loading
  Pajek files with vertex labels and edge weights.
*/

class PajekReader
{
public:

  template <class SimplicialComplex> void operator()( const std::string& filename, SimplicialComplex& K )
  {
    std::ifstream in( filename );
    if( !in )
      throw std::runtime_error( "Unable to read input file" );

    this->operator()( in, K );
  }

  template <class SimplicialComplex> void operator()( std::ifstream& in, SimplicialComplex& K )
  {
    using namespace aleph::utilities;

    using Simplex           = typename SimplicialComplex::ValueType;
    using DataType          = typename Simplex::DataType;
    using VertexType        = typename Simplex::VertexType;

    Mode mode = Mode::Unspecified;

    std::regex reKeyword  = std::regex( "\\*([[:alpha:]]+).*" );
    std::regex reKeyValue = std::regex( "\\*([[:alpha:]]+)[[:space:]]+([[:digit:]]+)" );
    std::regex reVertex   = std::regex( "([[:digit:]]+)[[:space:]]+\"(.*)\"");
    std::regex reEdge     = std::regex( "([[:digit:]]+)[[:space:]]+([[:digit:]]+)[[:space:]]*([[:digit:]\\.]+)*");

    std::smatch matches;

    std::vector<Simplex> vertices;
    std::vector<Simplex> edges;

    std::string line;
    while( std::getline( in, line ) )
    {
      // Although this is not explicitly specified in the rather terse
      // Pajek file format specification, it makes sense to remove all
      // white space characters in order to simplify parsing.
      line = trim( line );

      // Skip comments
      if( line.front() == '%' )
        continue;

      // 1st case: Found a keyword. This changes the parser mode and
      // requires us to load additional information.
      if( std::regex_match( line, matches, reKeyword ) )
      {
        auto name = matches[1];
        std::transform( name.begin(), name.end(), name.begin(), ::tolower );

        if( name == "vertices" || name == "verts" )
        {
          if( !std::regex_match( line, matches, reKeyValue ) )
            throw std::runtime_error( "Unable to parse vertices specification" );

          auto value = std::toul( matches[2] );
          mode       = Mode::Vertices;

          vertices.reserve( value );
        }
        else if( name == "edges" || name == "arcs" )
          mode = Mode::Edges;
      }

      // 2nd case: Proceed according to parser mode: vertices
      else if( mode == Mode::Vertices )
      {
        if( !std::regex_match( line, matches, reVertex ) )
          throw std::runtime_error( "Unable to parse vertex identifier" );

        auto id    = matches[1];
        auto label = matches[2];

        if( _labels.find(id) != _labels.end() )
          throw std::runtime_error( "Duplicate vertex identifier" );

        // Note that this construction is more generic than actually
        // specified in the description of the file format. This one
        // permits the use of non-contiguous node IDs.
        _labels[id] = label;

        // TODO: We could potentially support vertex weights here as well but
        // the original file format specification apparently does not account
        // for it.
        vertices.push_back( Simplex( VertexType( std::stoul(id) ) ) );
      }

      // 2nd case: Proceed according to parser mode: edges
      else if( mode == Mode::Edges )
      {
      }

      // 2nd case: Proceed according to parse mode: unspecified
      else if( mode == Mode::Unspecified )
      {
        // Happily ignore the line and churn along. It might make sense to
        // put up warnings here, though...
      }
    }
  }

private:

  enum class Mode
  {
    Vertices,
    Edges,
    Unspecified
  };

  // Local storage -----------------------------------------------------
  //
  // Variables in this section contain the result of the last parsing
  // process. This is useful when clients are querying attributes.

  std::map<std::string, std::string> _labels;
};

} // namespace io

} // namespace topology

} // namespace aleph

#endif
