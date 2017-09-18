#ifndef ALEPH_TOPOLOGY_IO_GML_HH__
#define ALEPH_TOPOLOGY_IO_GML_HH__

#include <fstream>
#include <map>
#include <set>
#include <regex>
#include <sstream>
#include <stack>
#include <string>
#include <unordered_map>
#include <vector>

#include <aleph/utilities/String.hh>

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

  - \c id (for nodes)
  - \c label (for nodes)
  - \c source (for edges)
  - \c target (for edges)
  - \c weight (for edges)
*/

class GMLReader
{
public:

  /**
    Reads a simplicial complex from a file, using the default maximum
    functor for weight assignment. If you want to change the functor,
    please refer to the overloaded variant of this method.

    @param filename Input filename
    @param  K       Simplicial complex

    @see operator()( const std::string&, SimplicialComplex& )
  */

  template <class SimplicialComplex> void operator()( const std::string& filename, SimplicialComplex& K )
  {
    using Simplex    = typename SimplicialComplex::ValueType;
    using DataType   = typename Simplex::DataType;

    this->operator()( filename, K, [] ( DataType a, DataType b ) { return std::max(a,b); } );
  }

  /**
    Reads a simplicial complex from a file while supporting arbitrary
    functors for weight assignment. The functor needs to support this
    interface:

    \code{.cpp}
    using SimplexType = typename SimplicialComplex::ValueType;
    using DataType    = typename Simplex::DataType;

    DataType Functor::operator()( DataType a, DataType b )
    {
      // Do something with a and b. Typically, this would be calculating
      // either the minimum or the maximum...
      return std::max(a, b);
    }
    \endcode

    Please refer to the documentation of SimplicialComplexReader::operator()( const std::string& SimplicialComplex&, Functor )
    for more details.
  */

  template <class SimplicialComplex, class Functor> void operator()( const std::string& filename, SimplicialComplex& K, Functor f )
  {
    std::ifstream in( filename );
    if( !in )
      throw std::runtime_error( "Unable to read input file" );

    this->operator()( in, K, f );
  }

  /** @overload operator()( const std::string&, SimplicialComplex& ) */
  template <class SimplicialComplex> void operator()( std::ifstream& in, SimplicialComplex& K )
  {
    using Simplex    = typename SimplicialComplex::ValueType;
    using DataType   = typename Simplex::DataType;

    this->operator()( in, K, [] ( DataType a, DataType b ) { return std::max(a,b); } );
  }

  /** @overload operator()( const std::string&, SimplicialComplex&, SimplicialComplex&, Functor ) */
  template <class SimplicialComplex, class Functor> void operator()( std::ifstream& in, SimplicialComplex& K, Functor f )
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
      line        = trim( line );
      auto tokens = split( line );

      // Skip comment lines. This is somewhat wasteful because I am
      // splitting the complete string even though I merely need to
      // know the first token.
      if( tokens.empty() == false && ( tokens.front() == "comment" || tokens.front() == "Creator" ) )
        continue;

      // A new level may also be opened "inline" by specifying a name
      // and the opening bracket at the same time.
      bool newLevel = isLevel( line ) || ( tokens.size() == 2 && tokens.back() == "[" && isLevel( tokens.front() ) );

      // Detecting a new level
      if( newLevel )
      {
        auto level = tokens.size() == 2 ? tokens.front() : line;

        if( lastLevel.empty() )
          lastLevel = level;
        else
          throw std::runtime_error( "Encountered incorrectly-nested levels" );

        // Only store the level directly when it has been specified
        // inline
        if( tokens.size() == 2 )
        {
          currentLevel.push( level );
          lastLevel = "";
        }
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
        if( currentLevel.empty() )
          throw std::runtime_error( "Expected a non-empty current level" );

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

    std::vector<Simplex> simplices;
    simplices.reserve( _nodes.size() + _edges.size() );

    // Maps an ID directly to its corresponding simplex in order to
    // facilitate the assignment of weights.
    std::unordered_map<VertexType, Simplex> id_to_simplex;

    for( auto&& node : _nodes )
    {
      auto id = static_cast<VertexType>( getID( nodeIDs, node.id ) );

      if( node.dict.find( "weight" ) != node.dict.end() )
        simplices.push_back( Simplex( id, convert<DataType>( node.dict.at( "weight" ) ) ) );
      else if( node.dict.find( "value" ) != node.dict.end() )
        simplices.push_back( Simplex( id, convert<DataType>( node.dict.at( "value" ) ) ) );
      else
        simplices.push_back( Simplex( id ) );

      id_to_simplex[id] = simplices.back();
    }

    // Create edges ----------------------------------------------------

    for( auto&& edge : _edges )
    {
      auto u = static_cast<VertexType>( getID( nodeIDs, edge.source ) );
      auto v = static_cast<VertexType>( getID( nodeIDs, edge.target ) );

      // No optional data attached; need to create weight based on node
      // weights, if those are available.
      if( edge.dict.find( "weight" ) == edge.dict.end() && edge.dict.find( "value" ) == edge.dict.end() )
      {
        auto uSimplex = id_to_simplex.at(u);
        auto vSimplex = id_to_simplex.at(v);

        simplices.push_back( Simplex( {u,v}, f( uSimplex.data(), vSimplex.data() ) ) );
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

  /**
    Returns a map that maps a node ID to an index. It corresponds
    to the order in which node IDs were allocated. It is *always*
    zero-indexed.
  */

  template <class VertexType> std::map<std::string, VertexType> id_to_index() const noexcept
  {
    std::set<std::string> nodeIDs;
    for( auto&& node : _nodes )
      nodeIDs.insert( node.id );

    std::map<std::string, VertexType> result;

    for( auto&& node : _nodes )
      result[ node.id ] = static_cast<VertexType>( getID( nodeIDs, node.id ) );

    return result;
  }

private:

  /**
    Auxiliary function for creating a numerical ID out of a parsed ID.
    In case non-numerical IDs are being used in the source file, this
    function ensures that internal IDs *always* start with a zero. If,
    however, numerical IDs are being used, they are converted as-is.
  */

  static std::size_t getID( const std::set<std::string>& ids, const std::string& id )
  {
    // Try to be smart: If the node ID can be converted into the
    // vertex type, we use the converted number instead.
    try
    {
      auto convertedID = std::stoll( id );
      return static_cast<std::size_t>( convertedID );
    }
    catch( std::out_of_range& e )
    {
      throw;
    }
    catch( std::invalid_argument& )
    {
      return static_cast<std::size_t>( std::distance( ids.begin(), ids.find( id ) ) );
    }
  };

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
