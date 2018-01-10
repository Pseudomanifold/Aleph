#ifndef ALEPH_TOPOLOGY_IO_GRAPHML_HH__
#define ALEPH_TOPOLOGY_IO_GRAPHML_HH__

namespace aleph
{

namespace topology
{

namespace io
{

/**
  @class GraphMLReader
  @brief Parses files in GraphML format

  This is a simple reader for graphs in GraphML format. Only a basic subset of
  the specification is supported, viz. reading nodes and edges, and extracting
  user-specified data.
*/

class GraphMLReader
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
  }

  /** Retrieves attribute names for the node attributes. */
  std::vector<std::string> getNodeAttributeNames() const
  {
    std::set<std::string> names;

    for( auto&& node : _nodes )
    {
      for( auto&& pair : node.dict )
        names.insert( pair.first );
    }

    return { names.begin(), names.end() };
  }

  /** Retrieves attribute names for the edge attributes */
  std::vector<std::string> getEdgeAttributeNames() const
  {
    std::set<std::string> names;

    for( auto&& edge : _edges )
    {
      for( auto&& pair : edge.dict )
        names.insert( pair.first );
    }

    return { names.begin(), names.end() };
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
      result[ node.id ] = static_cast<VertexType>( std::distance( nodeIDs.begin(), nodeIDs.find( node.id ) ) );

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
    bool isDirected;
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
    std::string source; // source node ID
    std::string target; // target node ID

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

} // namespace io

} // namespace topology

} // namespace aleph

#endif
