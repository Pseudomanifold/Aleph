#ifndef ALEPH_TOPOLOGY_IO_SPARSE_ADJACENCY_MATRIX_HH__
#define ALEPH_TOPOLOGY_IO_SPARSE_ADJACENCY_MATRIX_HH__

#include <aleph/topology/filtrations/Data.hh>

#include <aleph/utilities/Filesystem.hh>
#include <aleph/utilities/String.hh>

#include <algorithm>
#include <fstream>
#include <set>
#include <string>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <tuple>
#include <vector>

namespace aleph
{

namespace topology
{

namespace io
{

/**
  @class SparseAdjacencyMatrixReader
  @brief Parses files in sparse adjacency matrix format

  The purpose of this class is to read complex networks and represent
  them as simplicial complexes. The networks are assumed to be in the
  *sparse adjacency matrix* format developed by Dortmund University.

  The class has various optional settings for changing the way things
  are parsed and/or selecting different kinds of information that has
  to be retrieved from the graph.

  @see https://ls11-www.cs.tu-dortmund.de/staff/morris/graphkerneldatasets
  @see sparse_adjacency_matrices.cc
*/

class SparseAdjacencyMatrixReader
{
public:
  template <class SimplicialComplex> void operator()( const std::string& filename,
                                                      std::vector<SimplicialComplex>& complexes )
  {
    using Simplex    = typename SimplicialComplex::ValueType;
    using VertexType = typename Simplex::VertexType;
    using Edge       = std::pair<VertexType, VertexType>;

    std::unordered_set<VertexType> vertices;
    std::vector<Edge> edges;

    std::tie( vertices, edges )
      = this->readVerticesAndEdges<VertexType>( filename );

    auto graphIndicatorFilename = getFilenameGraphIndicator( filename );

    if( !aleph::utilities::exists( graphIndicatorFilename ) )
      throw std::runtime_error( "Missing required graph indicator file" );

    // Stores *all* graph IDs. The set makes sense here because I want
    // to ensure that repeated calls to this function always yield the
    // same order.
    std::set<VertexType> graphIDs;

    // Maps a node ID to a graph ID, i.e. yields the index of the graph
    // that should contain the node. All indices are starting at 1, but
    // will be mapped to 0 later on.
    std::unordered_map<VertexType, VertexType> node_id_to_graph_id;

    std::tie( graphIDs, node_id_to_graph_id )
      = readGraphAndNodeIDs<VertexType>( graphIndicatorFilename );

    using IndexType = std::size_t;

    // Maps a graph ID (arbitrary start point) to an index in the
    // vector.
    std::unordered_map<VertexType, IndexType> graph_id_to_index;

    {
      IndexType index = IndexType();

      for( auto&& id : graphIDs )
        graph_id_to_index[id] = index++;
    }

    // Maps a node ID (arbitrary start point) to an index. Note that no
    // direct vector for nodes/vertices exists; this is merely required
    // for internal bookkeeping.
    //
    // Moreover, this information can be used to access node attributes
    // directly.

    std::unordered_map<VertexType, IndexType> node_id_to_index;

    {
      std::set<VertexType> vertices_( vertices.begin(), vertices.end() );

      // Add also those nodes that are only defined implicitly because
      // they are isolated. This requires traversing the map that maps
      // node IDs to graph IDs.
      for( auto&& pair : node_id_to_graph_id )
        vertices_.insert( pair.first );

      IndexType index = IndexType();

      for( auto&& id : vertices_ )
        node_id_to_index[id] = index++;
    }

    // Reading optional attributes -------------------------------------

    if( _readGraphLabels )
      this->readGraphLabels( filename );

    if( _readNodeLabels )
      this->readNodeLabels( filename );

    if( _readNodeAttributes )
      this->readNodeAttributes( filename );

    if( _readEdgeAttributes )
      this->readEdgeAttributes( filename );

    // Create output ---------------------------------------------------
    //
    // Create the set of output graphs and distribute the edges among
    // them according to their graph ID. This function also does some
    // sanity checks in order to check input data consistency.

    complexes.clear();
    complexes.resize( graphIDs.size() );

    // Contains the *actual* labels of all graphs. It is possible that
    // some of the input data files do not contain a contiguous series
    // of labels, making it necessary to store the ones that have been
    // encountered.
    std::vector<std::string> labels( graphIDs.size() );

    using DataType = typename Simplex::DataType;

    for( auto&& vertex : vertices )
    {
      auto&& id    = node_id_to_graph_id[vertex];
      auto&& index = graph_id_to_index[id];
      auto&& K     = complexes[index];
      auto s       = Simplex( vertex );

      // Most of the time, this just degenerates into an identity
      // lookup, but in case labels are not contiguous, they will
      // be assigned correctly here.
      labels[index] = _graphLabels.at( index );

      if( _readNodeAttributes && isValidIndex( _nodeAttributeIndex ) )
      {
        auto&& index = node_id_to_index[vertex];
        s.setData( static_cast<DataType>( _nodeAttributes[index][ _nodeAttributeIndex ] ) );
      }

      K.push_back( s );
    }

    for( std::size_t i = 0; i < edges.size(); i++ )
    {
      auto&& edge = edges.at(i);
      auto&& u    = edge.first;
      auto&& v    = edge.second;
      auto&& uID  = node_id_to_graph_id[u];
      auto&& vID  = node_id_to_graph_id[v];

      if( uID != vID )
        throw std::runtime_error( "Format error: an edge must not belong to multiple graphs" );

      auto&& index = graph_id_to_index[ uID ];
      auto&& K     = complexes[index];
      auto s       = Simplex( {u,v} );

      if( _readEdgeAttributes && isValidIndex( _edgeAttributeIndex ) )
        s.setData( static_cast<DataType>( _edgeAttributes[i][ _edgeAttributeIndex ] ) );

      K.push_back( s );
    }

    for( auto&& K : complexes )
      K.sort( aleph::topology::filtrations::Data<Simplex>() );

    if( labels.size() < _graphLabels.size() )
      _graphLabels = labels;
  }

  // Output ------------------------------------------------------------

  /**
    Stores graph labels in an output iterator. This function does not
    check for their availability.
  */

  template <class OutputIterator> void graphLabels( OutputIterator result )
  {
    std::copy( _graphLabels.begin(), _graphLabels.end(), result );
  }

  /**
    Stores node labels in an output iterator. This function does not
    check for their availability.
  */

  template <class OutputIterator> void nodeLabels( OutputIterator result )
  {
    std::copy( _graphLabels.begin(), _graphLabels.end(), result );
  }

  // Configuration options ---------------------------------------------
  //
  // The following attributes configure how the parsing process works
  // and which attributes are being read.

  void setSeparator( const std::string& separator ) noexcept
  {
    // I am not performing any sanity checks here. If the separator does
    // not constitute a valid regular expression, the parsing process is
    // just going to fail later on.
    _separator = separator;
  }

  std::string separator() const noexcept
  {
    return _separator;
  }

  void setReadGraphLabels( bool value = true )    noexcept { _readGraphLabels = value;    }
  void setReadNodeLabels( bool value = true )     noexcept { _readNodeLabels = value;     }
  void setReadNodeAttributes( bool value = true ) noexcept { _readNodeAttributes = value; }
  void setReadEdgeAttributes( bool value = true ) noexcept { _readEdgeAttributes = value; }
  void setTrimLines( bool value = true )          noexcept { _trimLines = value;          }

  void setNodeAttributeIndex( std::size_t value ) noexcept { _nodeAttributeIndex = value; }
  void setEdgeAttributeIndex( std::size_t value ) noexcept { _edgeAttributeIndex = value; }

  bool readGraphLabels()    const noexcept { return _readGraphLabels;    }
  bool readNodeLabels()     const noexcept { return _readNodeLabels;     }
  bool readNodeAttributes() const noexcept { return _readNodeAttributes; }
  bool readEdgeAttributes() const noexcept { return _readEdgeAttributes; }
  bool trimLines()          const noexcept { return _trimLines;          }

  std::size_t nodeAttributeIndex() const noexcept { return _nodeAttributeIndex; }
  std::size_t edgeAttributeIndex() const noexcept { return _edgeAttributeIndex; }

private:

  /**
    Reads all vertices and  edges from a sparse adjacency matrix. Note
    that this function is not static because it needs access to values
    that depend on class instances.
  */

  template <class VertexType>
    std::pair<
      std::unordered_set<VertexType>,
      std::vector< std::pair<VertexType, VertexType> >
    > readVerticesAndEdges( const std::string& filename )
  {
    using Edge = std::pair<VertexType, VertexType>;

    std::unordered_set<VertexType> vertices;
    std::vector<Edge> edges;

    std::ifstream in( filename );
    if( !in )
      throw std::runtime_error( "Unable to read input adjacency matrix file" );

    std::string line;

    while( std::getline( in, line ) )
    {
      using namespace aleph::utilities;

      auto tokens = split( line, _separator );

      if( tokens.size() == 2 )
      {
        auto u = convert<VertexType>( tokens.front() );
        auto v = convert<VertexType>( tokens.back()  );

        edges.push_back( std::make_pair(u, v) );

        vertices.insert( u );
        vertices.insert( v );
      }
      else
        throw std::runtime_error( "Format error: cannot parse line in sparse adjacency matrix" );
    }

    return std::make_pair( vertices, edges );
  }

  /**
    Reads graph and node IDs from an input file. The node ID is
    implicitly encoded by the current line number. Each line in
    turn contains a graph identifier. This is usually a number,
    but the function does not require IDs to be contiguous.

    @param filename Input filename

    @returns Pair that contains information stored in the graph. The
    first entry contains the set of *all* graph IDs. The second item
    contains a dictionary for mapping a node ID to its corresponding
    graph ID.
  */

  template <class VertexType>
    std::pair<
      std::set<VertexType>,
      std::unordered_map<VertexType, VertexType>
    > readGraphAndNodeIDs( const std::string& filename ) const
  {
    std::unordered_map<VertexType, VertexType> node_id_to_graph_id;
    std::set<VertexType> graphIDs;

    std::ifstream in( filename );
    if( !in )
      throw std::runtime_error( "Unable to read graph indicator file" );

    std::string line;
    VertexType nodeID = static_cast<VertexType>( _firstNodeID );

    while( std::getline( in, line ) )
    {
      using namespace aleph::utilities;
      line = trim( line );

      bool success                    = false;
      auto graphID                    = convert<VertexType>( line, success );

      if( !success )
        throw std::runtime_error( "Unable to convert graph ID to numerical type" );

      node_id_to_graph_id[ nodeID++ ] = graphID;

      graphIDs.insert( graphID );
    }

    return std::make_pair( graphIDs, node_id_to_graph_id );
  }

  std::vector<std::string> readLabels( const std::string& filename )
  {
    std::ifstream in( filename );
    if( !in )
      throw std::runtime_error( "Unable to read labels input file" );

    std::vector<std::string> labels;
    std::string line;

    while( std::getline( in, line ) )
    {
      if( _trimLines )
        line = aleph::utilities::trim( line );

      labels.push_back( line );
    }

    return labels;
  }

  void readGraphLabels( const std::string& filename )
  {
    auto graphLabelsFilename = getFilenameGraphLabels( filename );
    _graphLabels             = readLabels( graphLabelsFilename );
  }

  void readNodeLabels( const std::string& filename )
  {
    auto nodeLabelsFilename = getFilenameNodeLabels( filename );
    _nodeLabels             = readLabels( nodeLabelsFilename );
  }

  std::vector< std::vector<double> > readAttributes( const std::string& filename )
  {
    std::ifstream in( filename );
    if( !in )
      throw std::runtime_error( "Unable to read attributes input file" );

    std::vector< std::vector<double> > allAttributes;

    std::string line;
    while( std::getline( in, line ) )
    {
      using namespace aleph::utilities;

      if( _trimLines )
        line = trim( line );

      auto tokens = split( line, _separator );

      std::vector<double> attributes;
      std::transform( tokens.begin(), tokens.end(), std::back_inserter( attributes ),
        [] ( const std::string& token )
        {
          return convert<double>( token );
        }
      );

      attributes.shrink_to_fit();
      allAttributes.push_back( attributes );
    }

    return allAttributes;
  }

  void readNodeAttributes( const std::string& filename )
  {
    auto nodeAttributesFilename = getFilenameNodeAttributes( filename );
    _nodeAttributes             = readAttributes( nodeAttributesFilename );
  }

  void readEdgeAttributes( const std::string& filename )
  {
    auto edgeAttributesFilename = getFilenameEdgeAttributes( filename );
    _edgeAttributes             = readAttributes( edgeAttributesFilename );
  }

  static bool isValidIndex( std::size_t index )
  {
    return index != std::numeric_limits<std::size_t>::max();
  }

  /**
   Given a base filename, gets its prefix. The prefix is everything that
   comes before the last `_` character. It is used to generate filenames
   for various auxiliary information about the graph if not specified by
   the client.
  */

  static std::string getPrefix( const std::string& filename )
  {
    // Note that this contains the complete filename portion of the path
    // along with any subdirectories.
    auto prefix = filename.substr( 0, filename.find_last_of( '_' ) );
    return prefix;
  }

  /*
    Given a base filename, calculates the default filename for the graph
    indicator values, for example. This filename is generated by getting
    a prefix of the filename and attaching e.g. `_graph_indicator.txt`.

    This function is *not* used if the user specified a file manually in
    which to find such information.
  */

  static std::string getFilenameGraphIndicator( const std::string& filename )  { return getPrefix(filename) + "_graph_indicator.txt";  }
  static std::string getFilenameGraphLabels( const std::string& filename )     { return getPrefix(filename) + "_graph_labels.txt";     }
  static std::string getFilenameNodeLabels( const std::string& filename )      { return getPrefix(filename) + "_node_labels.txt";      }
  static std::string getFilenameEdgeLabels( const std::string& filename )      { return getPrefix(filename) + "_edge_labels.txt";      }
  static std::string getFilenameEdgeAttributes( const std::string& filename )  { return getPrefix(filename) + "_edge_attributes.txt";  }
  static std::string getFilenameNodeAttributes( const std::string& filename )  { return getPrefix(filename) + "_node_attributes.txt";  }
  static std::string getFilenameGraphAttributes( const std::string& filename ) { return getPrefix(filename) + "_graph_attributes.txt"; }

  bool _readGraphLabels    = true;
  bool _readNodeLabels     = false;
  bool _readNodeAttributes = false;
  bool _readEdgeAttributes = false;
  bool _trimLines          = true;

  std::size_t _nodeAttributeIndex = std::numeric_limits<std::size_t>::max();
  std::size_t _edgeAttributeIndex = std::numeric_limits<std::size_t>::max();

  /**
    By default, the first node ID starts with 1. This should be
    changeable, however, because we may just as well have files
    for which the initial node ID is zero-based.
  */

  std::size_t _firstNodeID = 1;

  /**
    Graph labels stored during the main parsing routine of this class.
    If no graph labels are specified, this vector remains empty. Since
    the format does not specify the format of graph labels, the labels
    will *not* be converted but reported as-is.

    Graph labels are only read if `_readGraphLabels` is true.
  */

  std::vector<std::string> _graphLabels;

  /** Node labels; the same comments as for the graph labels apply */
  std::vector<std::string> _nodeLabels;

  /**
    Node attributes stored during the main parsing routine of this
    class. It is assumed that the attributes are *numerical* or at
    least *convertible* to a numerical representation.
  */

  std::vector< std::vector<double> > _nodeAttributes;

  /**
    Edge attributes; just like the node attributes, it is assumed that
    they are convertible to a numerical representation.
  */

  std::vector< std::vector<double> > _edgeAttributes;

  /** Default separator between edges. This works strings such as `1,2`. */
  std::string _separator = ",";
};

} // namespace io

} // namespace topology

} // namespace aleph

#endif
