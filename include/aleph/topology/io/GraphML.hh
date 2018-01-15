#ifndef ALEPH_TOPOLOGY_IO_GRAPHML_HH__
#define ALEPH_TOPOLOGY_IO_GRAPHML_HH__

#include <aleph/config/TinyXML2.hh>

#include <aleph/utilities/String.hh>

#ifdef ALEPH_WITH_TINYXML2
  #include <tinyxml2.h>
#endif

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
    Default dictionary type used by this class. This is a convenience
    declaration that shortens the code.
  */

  using Dictionary = std::map<std::string, std::string>;

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
    _graph = {};

    _nodes.clear();
    _edges.clear();

    #ifdef ALEPH_WITH_TINYXML2
      using namespace tinyxml2;

      XMLDocument document;
      document.LoadFile( filename.c_str() );

      auto root        = document.FirstChild();
      auto graphml     = root->NextSiblingElement( "graphml" );

      // 1. Read optional information about keys in the graph ----------

      {
        auto key = graphml->FirstChildElement( "key" );
        while( key )
        {
          auto id   = key->Attribute( "id" );
          auto name = key->Attribute( "attr.name" );
          auto type = key->Attribute( "for" );

          // It makes sense to switch these two items because we use the
          // name, e.g. 'weight', to look up the ID, e.g. 'id0'.
          if( std::string( type ) == "node" )
            _graph.nodeKeys[name] = id;
          else if( std::string( type ) == "edge" )
            _graph.edgeKeys[name] = id;
          else
            throw std::runtime_error( "Attribute must belong to either nodes or edges" );

          key = key->NextSiblingElement( "key" );
        }
      }

      // 2. Start parsing the graph ------------------------------------

      auto graph = graphml->FirstChildElement( "graph" );
      if( !graph )
        throw std::runtime_error( "GraphML file has to contain at least one graph" );

      // Check whether the graph is directed or not. Note that this does
      // not have any influence on the parsing process for now.
      auto edgedefault = graph->Attribute( "edgedefault" );
      if( edgedefault && std::string( edgedefault ) == "directed" )
        _graph.isDirected = true;
      else
        _graph.isDirected = false;

      auto current = graph->FirstChildElement();

      while( current )
      {
        std::string name = current->Name();
        if( name == "node" )
          parseNode( current );
        else if( name == "edge" )
          parseEdge( current );
        else
        {
          // Ignoring *unknown* elements for now, because they could not
          // possibly harm the extraction process.
        }

        current = current->NextSiblingElement();
      }
    #endif

    // 3. Create simplicial complex ----------------------------------
    //
    // This works regardless of the presence of an XML parsing library,
    // even though the resulting complex may of course be empty.

    using Simplex    = typename SimplicialComplex::ValueType;
    using DataType   = typename Simplex::DataType;
    using VertexType = typename Simplex::VertexType;

    std::vector<Simplex> simplices;
    simplices.reserve( _nodes.size() + _edges.size() );

    auto id_to_index_map = id_to_index<VertexType>();

    using namespace aleph::utilities;

    for( auto&& node : _nodes )
    {
      auto vertex = id_to_index_map.at( node.id );
      auto weight = DataType();

      if( _readNodeWeights && not _nodeWeightAttribute.empty() )
      {
        if( _graph.nodeKeys.find( _nodeWeightAttribute ) != _graph.nodeKeys.end() )
        {
          bool success = false;
          auto data    = node.dict.at( _graph.nodeKeys.at( _nodeWeightAttribute ) );
          weight       = convert<DataType>( data, success );

          if( !success )
            throw std::runtime_error( "Unable to convert node weight to data type" );
        }
      }

      simplices.push_back( Simplex( vertex, weight ) );
    }
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

  void setReadNodeWeights( bool value = true ) noexcept           { _readNodeWeights = value; }
  void setReadEdgeWeights( bool value = true ) noexcept           { _readEdgeWeights = value; }

  bool readNodeWeights() const noexcept                           { return _readNodeWeights; }
  bool readEdgeWeights() const noexcept                           { return _readEdgeWeights; }

  void setNodeWeightAttribute( const std::string& name ) noexcept { _nodeWeightAttribute = name; }
  void setEdgeWeightAttribute( const std::string& name ) noexcept { _edgeWeightAttribute = name; }

  const std::string& nodeWeightAttribute() const noexcept         { return _nodeWeightAttribute; }
  const std::string& edgeWeightAttribute() const noexcept         { return _edgeWeightAttribute; }

private:

  #ifdef ALEPH_WITH_TINYXML2

    void parseNode( tinyxml2::XMLElement* element )
    {
      auto name = element->Name();
      auto id   = element->Attribute( "id" );

      if( std::string( name ) != "node" )
        throw std::runtime_error( "Unexpected element for node parsing" );

      if( !id )
        throw std::runtime_error( "Node element must specify ID" );

      Node node;
      node.id = id;

      // Parse additional details of the node, even though they may not
      // be used in the creation of a simplicial complex.
      parseData( element, node.dict );

      _nodes.push_back( node );
    }

    void parseEdge( tinyxml2::XMLElement* element )
    {
      auto name   = element->Name();
      auto id     = element->Attribute( "id" );
      auto source = element->Attribute( "source" );
      auto target = element->Attribute( "target" );

      if( std::string( name ) != "edge" )
        throw std::runtime_error( "Unexpected element for edge parsing" );

      if( !id )
        throw std::runtime_error( "Edge element must specify ID" );

      if( !source || !target )
        throw std::runtime_error( "Edge element must specify both source and target" );

      Edge edge;
      edge.source = source;
      edge.target = target;

      // Parse additional details of the edge, even though they may not
      // be used in the creation of a simplicial complex.
      parseData( element, edge.dict );
    }

    void parseData( tinyxml2::XMLElement* element, Dictionary& dict )
    {
      auto child = element->FirstChildElement();
      while( child )
      {
        std::string name = child->Name();
        if( name == "data" )
        {
          auto key   = child->Attribute( "key" );
          auto value = child->GetText();

          dict[key] = value;
        }

        child = child->NextSiblingElement();
      }
    }

  #endif

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

    // This maps all key IDs to their corresponding names, which is how
    // they are stored for nodes and edges.
    Dictionary nodeKeys;
    Dictionary edgeKeys;
  };

  /** Describes a parsed node along with all of its attributes */
  struct Node
  {
    std::string id;

    Dictionary dict; // all remaining attributes
  };

  /** Describes a parsed edge along with all of its attributes */
  struct Edge
  {
    std::string source; // source node ID
    std::string target; // target node ID

    Dictionary dict; // all remaining attributes
  };

  // Attributes --------------------------------------------------------

  bool _readEdgeWeights = true;
  bool _readNodeWeights = true;

  std::string _edgeWeightAttribute = "weight"; // default attribute for node weight extraction
  std::string _nodeWeightAttribute = "weight"; // default attribute for edge weight extraction

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
