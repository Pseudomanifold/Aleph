#ifndef ALEPH_TOPOLOGY_IO_GML_HH__
#define ALEPH_TOPOLOGY_IO_GML_HH__

#include <fstream>
#include <set>
#include <stack>
#include <string>
#include <vector>

// FIXME: Remove after debugging
#include <iostream>

#include "utilities/String.hh"

namespace aleph
{

namespace topology
{

namespace io
{

/**
  @class GMLReader
  @brief Parses files in GML (Graph Modeling Language) format

  This is a simple reader for graphs in GML format. It suports a basic
  subset of the GML specification, viz. the specification of different
  attributes for nodes, as well as weight specifications for edges.

  Currently, the following attributes will be read:

  * \c id (for nodes)
  * \c label (for nodes)
  * \c source (for edges)
  * \c target (for edges)
  * \c weight (for edges)
*/

class GMLReader
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

    std::string line;
    std::set<std::string> levels = { "graph", "node", "edge" };

    auto isLevel = [&levels] ( const std::string& name )
    {
      return levels.find( name ) != levels.end();
    };

    // Specifies the current level the parse is in. May be either one of
    // the known levels above.
    std::stack<std::string> currentLevel;

    // Last level that was read by the parser. If an open bracket '[' is
    // identified, this will become the current level.
    std::string lastLevel;

    while( std::getline( in, line ) )
    {
      line = trim( line );


      // Skip comment lines. This is somewhat wasteful because I am
      // splitting the complete string even though I merely need to
      // know the first token.
      {
        auto tokens = split( line );
        if( tokens.empty() == false && tokens.front() == "comment" )
          continue;
      }

      // Detecting a new level
      if( isLevel( line ) )
      {
        if( lastLevel.empty() )
          lastLevel = line;
        else
          throw std::runtime_error( "Encountered incorrectly-nested levels" );
      }

      // Opening a new level
      else if( line == "[" )
      {
        std::cerr << "* Entering level = " << lastLevel << "\n";
        currentLevel.push( lastLevel );
        lastLevel = "";
      }

      // Closing a new level
      else if( line == "]" )
      {
        std::cerr << "* Leaving level = " << currentLevel.top() << "\n";
        currentLevel.pop();
      }
    }

    (void) K;
  }
};

} // namespace io

} // namespace topology

} // namespace aleph

// graph
// [
//   node
//   [
//    id A
//   ]
//   node
//   [
//    id B
//   ]
//   node
//   [
//    id C
//   ]
//    edge
//   [
//    source B
//    target A
//   ]
//   edge
//   [
//    source C
//    target A
//   ]
// ]

#endif
