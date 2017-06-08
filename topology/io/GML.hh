#ifndef ALEPH_TOPOLOGY_IO_GML_HH__
#define ALEPH_TOPOLOGY_IO_GML_HH__

#include <fstream>
#include <map>
#include <set>
#include <regex>
#include <sstream>
#include <stack>
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
    _nodes.clear();
    _edges.clear();

    using namespace aleph::utilities;

    using Simplex           = typename SimplicialComplex::ValueType;
    using DataType          = typename Simplex::DataType;
    using VertexType        = typename Simplex::VertexType;

    std::string line;
    std::set<std::string> levels     = { "graph", "node", "edge" };
    std::set<std::string> attributes = { "id", "label", "source", "target", "value", "weight" };

    auto isLevel = [&levels] ( const std::string& name )
    {
      return levels.find( name ) != levels.end();
    };

    auto isAttribute = [&attributes] ( const std::string& name )
    {
      return attributes.find( name ) != attributes.end();
    };

    // Specifies the current level the parser is in. May be either one
    // of the known levels above.
    std::stack<std::string> currentLevel;

    // Last level that was read by the parser. If an open bracket '[' is
    // identified, this will become the current level.
    std::string lastLevel;

    Graph graph;
    Node node;
    Edge edge;

    std::regex reAttribute = std::regex( "([[:alpha:]]+)[[:space:]]*.*" );
    std::regex reKeyValue  = std::regex( "([[:alpha:]]+)[[:space:]]+([[:alnum:]\\.]+)" );
    std::regex reLabel     = std::regex( "(label)[[:space:]]+\"([^\"]+)\"" );

    while( std::getline( in, line ) )
    {
      line = trim( line );

      // Skip comment lines. This is somewhat wasteful because I am
      // splitting the complete string even though I merely need to
      // know the first token.
      {
        auto tokens = split( line );
        if( tokens.empty() == false && ( tokens.front() == "comment" || tokens.front() == "Creator" ) )
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
        currentLevel.push( lastLevel );
        lastLevel = "";
      }

      // Closing a new level
      else if( line == "]" )
      {
        if( currentLevel.top() == "node" )
          _nodes.push_back( node );
        else if( currentLevel.top() == "edge" )
          _edges.push_back( edge );

        // Reset node and edge data structure to fill them again once
        // a new level is being encountered.
        node = {};
        edge = {};

        currentLevel.pop();
      }

      // Check for attributes
      else
      {
        std::smatch matches;

        auto* dict = currentLevel.top() == "node" ? &node.dict
                                                  : currentLevel.top() == "edge" ? &edge.dict
                                                                                 : currentLevel.top() == "graph" ? &graph.dict
                                                                                                                 : throw std::runtime_error( "Current level is unknown" );

        if( std::regex_match( line, matches, reAttribute ) )
        {
          auto name = matches[1];
          if( isAttribute( name ) )
          {
            // Special matching for labels
            if( name == "label" )
              std::regex_match( line, matches, reLabel );

            // Regular matching for all other attributes
            else
              std::regex_match( line, matches, reKeyValue );

            auto value = matches[2];

            if( name == "id" )
              node.id = value;
            else if( name == "source" )
              edge.source = value;
            else if( name == "target" )
              edge.target = value;

            // Just add it to the dictionary of optional values
            else
             dict->operator[]( name ) = value;
          }
          // Skip unknown attributes...
          else
          {
          }
        }
      }
    }

    // Creates nodes (vertices) ----------------------------------------

    std::set<std::string> nodeIDs;

    for( auto&& node : _nodes )
    {
      auto pair = nodeIDs.insert( node.id );
      if( !pair.second )
        throw std::runtime_error( "Duplicate node id '" + node.id + "'" );
    }

    // Lambda expression for creating a numerical ID out of a parsed ID. In
    // case non-numerical IDs are being used in the source file this Lambda
    // expression ensures that internal IDs always start with a zero.
    auto getID = [&nodeIDs] ( const std::string& id )
    {
      // Try to be smart: If the node ID can be converted into the
      // vertex type, we use the converted number instead.
      try
      {
        auto convertedID = std::stoll( id );
        return static_cast<VertexType>( convertedID );
      }
      catch( std::out_of_range& e )
      {
        throw;
      }
      catch( std::invalid_argument& )
      {
        return static_cast<VertexType>(
          std::distance( nodeIDs.begin(), nodeIDs.find( id ) )
        );
      }
    };

    std::vector<Simplex> simplices;
    simplices.reserve( _nodes.size() + _edges.size() );

    for( auto&& node : _nodes )
    {
      auto id = getID( node.id );

      if( node.dict.find( "weight" ) != node.dict.end() )
        simplices.push_back( Simplex( id, convert<DataType>( node.dict.at( "weight" ) ) ) );
      else if( node.dict.find( "value" ) != node.dict.end() )
        simplices.push_back( Simplex( id, convert<DataType>( node.dict.at( "value" ) ) ) );
      else
        simplices.push_back( Simplex( id ) );
    }

    // Create edges ----------------------------------------------------

    auto getSimplexByID = [&simplices] ( VertexType id )
    {
      auto position = std::find_if( simplices.begin(), simplices.end(),
                                    [&id] ( const Simplex& s )
                                    {
                                      return s.dimension() == 0 && s[0] == id;
                                    } );

      if( position != simplices.end() )
        return *position;
      else
        throw std::runtime_error( "Querying unknown simplex for edge creation" );
    };

    for( auto&& edge : _edges )
    {
      auto u = getID( edge.source );
      auto v = getID( edge.target );

      // No optional data attached; need to create weight based on node
      // weights, if those are available.
      if( edge.dict.find( "weight" ) == edge.dict.end() && edge.dict.find( "value" ) == edge.dict.end() )
      {
        // TODO: This is very slow and could be sped up by maintaining
        // the simplex ID somewhere else with faster access.
        auto uSimplex = getSimplexByID( u );
        auto vSimplex = getSimplexByID( v );

        // TODO: Permit the usage of other weight assignment strategies
        // here, for example by using a functor.
        auto data = std::max( uSimplex.data(), vSimplex.data() );

        simplices.push_back( Simplex( {u,v}, data ) );
      }

      // Use converted weight
      else if( edge.dict.find( "weight" ) != edge.dict.end() )
        simplices.push_back( Simplex( {u,v}, convert<DataType>( edge.dict.at( "weight" ) ) ) );
      else if( edge.dict.find( "value" ) != edge.dict.end() )
        simplices.push_back( Simplex( {u,v}, convert<DataType>( edge.dict.at( "value" ) ) ) );
    }

    _graph = graph;
    K      = SimplicialComplex( simplices.begin(), simplices.end() );
  }

  /**
    Retrieves a map of attribute values for each node. The attributes
    will be assigned to the node ID. Empty attributes signify that no
    such information is available.
  */

  std::map<std::string, std::string> getNodeAttribute( const std::string& attribute ) const
  {
    std::map<std::string, std::string> map;

    for( auto&& node : _nodes )
    {
      if( node.dict.find( attribute ) == node.dict.end() )
        map[ node.id ] = std::string();
      else
        map[ node.id ] = node.dict.at( attribute );
    }

    return map;
  }

private:

  /** Describes a parsed graph along with all of its attributes */
  struct Graph
  {
    std::map<std::string, std::string> dict; // all remaining attributes
  };

  /** Describes a parsed node along with all of its attributes */
  struct Node
  {
    std::string id;

    std::map<std::string, std::string> dict; // all remaining attributes
  };

  /** Describes a parsed edge along with all of its attributes */
  struct Edge
  {
    std::string source;
    std::string target;

    std::map<std::string, std::string> dict; // all remaining attributes
  };

  // Local storage -----------------------------------------------------
  //
  // Variables in this section contain the result of the last parsing
  // process. This is useful when clients are querying attributes.

  Graph _graph;

  std::vector<Node> _nodes;
  std::vector<Edge> _edges;
};

/**
  @class GMLWriter
  @brief Writes files in GML (Graph Modeling Language) format

  This is a simple writer for graphs in GML format. It suports a basic
  subset of the GML specification, viz. the specification of different
  attributes for nodes, as well as weight specifications for edges.

  Given a simplicial complex, it will store it as a weighted graph.

  Currently, the following attributes will be written:

  * \c id (for nodes)
  * \c source (for edges)
  * \c target (for edges)
  * \c weight (for edges)
*/

class GMLWriter
{
public:

  template <class SimplicialComplex> void operator()( const std::string& filename, SimplicialComplex& K )
  {
    std::ofstream out( filename );
    if( !out )
      throw std::runtime_error( "Unable to open output file" );

    this->operator()( out, K );
  }

  template <class SimplicialComplex> void operator()( std::ostream& out, SimplicialComplex& K )
  {
    std::ostringstream streamNodes;
    std::ostringstream streamEdges;

    for( auto&& simplex : K )
    {
      if( simplex.dimension() == 0 )
      {
        streamNodes << "  node [\n"
                    << "    id " << *simplex.begin() << "\n"
                    << "  ]\n";
      }
      else if( simplex.dimension() == 1 )
      {
        auto u = *( simplex.begin()     );
        auto v = *( simplex.begin() + 1 );

        streamEdges << "  edge [\n"
                    << "    source " << u << "\n"
                    << "    target " << v << "\n"
                    << "    weight " << simplex.data() << "\n"
                    << "  ]\n";
      }
    }

    out << "graph [\n"
        << "  directed 0\n"
        << streamNodes.str() << "\n"
        << streamEdges.str() << "\n"
        << "]\n";
  }
};

} // namespace io

} // namespace topology

} // namespace aleph

#endif
