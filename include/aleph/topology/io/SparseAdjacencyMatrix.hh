#ifndef ALEPH_TOPOLOGY_IO_SPARSE_ADJACENCY_MATRIX_HH__
#define ALEPH_TOPOLOGY_IO_SPARSE_ADJACENCY_MATRIX_HH__

#include <aleph/utilities/String.hh>

// FIXME: remove after debugging
#include <iostream>

#include <fstream>
#include <set>
#include <string>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace aleph
{

namespace topology
{

namespace io
{

class SparseAdjacencyMatrixReader
{
public:
  template <class SimplicialComplex> void operator()( const std::string& filenameMatrix,
                                                      const std::string& filenameIndicator,
                                                      std::vector<SimplicialComplex>& result )
  {
    using Simplex    = typename SimplicialComplex::ValueType;
    using VertexType = typename Simplex::VertexType;
    using Edge       = std::pair<VertexType, VertexType>;

    std::vector<Edge> edges;
    std::unordered_set<VertexType> vertices;

    {
      std::ifstream in( filenameMatrix );
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
        {
          // TODO: throw error?
        }
      }

      std::cerr << __PRETTY_FUNCTION__ << ": Read " << edges.size() << " edges\n";
    }

    // Maps a node ID to a graph ID, i.e. yields the inex of the graph
    // that should contain the node. All indices are starting at 1 but
    // will be mapped to 0 later on.
    std::unordered_map<VertexType, VertexType> node_id_to_graph_id;

    // Stores *all* graph IDs. The set makes sense here because I want
    // to ensure that repeated calls to this function always yield the
    // same order.
    std::set<VertexType> graphIDs;

    {
      std::ifstream in( filenameIndicator );
      if( !in )
        throw std::runtime_error( "Unable to read input graph indicator file" );

      std::string line;
      VertexType nodeID = 1;

      while( std::getline( in, line ) )
      {
        using namespace aleph::utilities;
        line = trim( line );

        auto graphID                    = convert<VertexType>( line );
        node_id_to_graph_id[ nodeID++ ] = graphID;

        graphIDs.insert( graphID );
      }
    }

    // Maps a graph ID (arbitrary start point) to an index in the
    // vector.

    using IndexType = std::size_t;
    std::unordered_map<VertexType, IndexType> graph_id_to_index;

    {
      IndexType index = IndexType();

      for( auto&& id : graphIDs )
        graph_id_to_index[id] = index++;
    }

    // Create output ---------------------------------------------------
    //
    // Create the set of output graphs and distribute the edges among
    // them according to their graph ID. This function also does some
    // sanity checks in order to check input data consistency.

    result.clear();
    result.resize( graphIDs.size() );

    for( auto&& vertex : vertices )
    {
      auto&& id    = node_id_to_graph_id[vertex];
      auto&& index = graph_id_to_index[id];
      auto&& K     = result[index];

      K.push_back( Simplex( vertex ) );
    }

    for( auto&& edge : edges )
    {
      auto&& u   = edge.first;
      auto&& v   = edge.second;
      auto&& uID = node_id_to_graph_id[u];
      auto&& vID = node_id_to_graph_id[v];

      if( uID != vID )
        throw std::runtime_error( "Format error: an edge must not belong to multiple graphs" );

      auto&& index = graph_id_to_index[ uID ];
      auto&& K     = result[index];

      K.push_back( Simplex( {u,v} ) );
    }
  }

private:

  // TODO: make configurable
  std::string _separator = ",";
};

} // namespace io

} // namespace topology

} // namespace aleph

#endif
