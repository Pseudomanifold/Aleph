#ifndef ALEPH_PERSISTENT_HOMOLOGY_EXTENDED_PERSISTENCE_HIERARCHY__
#define ALEPH_PERSISTENT_HOMOLOGY_EXTENDED_PERSISTENCE_HIERARCHY__

#include <boost/bimap.hpp>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/breadth_first_search.hpp>

#include "topology/SimplicialComplex.hh"

namespace aleph
{

namespace detail
{

using AdjacencyGraph = boost::adjacency_list<
  boost::vecS,
  boost::vecS,
  boost::undirectedS,
  boost::no_property,
  boost::no_property>;

/** Extracts the adjacency graph of a simplicial complex */
template <class Simplex>
std::pair<boost::bimap<typename Simplex::vertex_type, unsigned>, AdjacencyGraph> extractZeroDimensionalAdjacencyGraph( const topology::SimplicialComplex<Simplex>& S )
{
  AdjacencyGraph adjacencyGraph;

  using Vertex           = typename Simplex::vertex_type;
  using VertexDescriptor = boost::graph_traits<AdjacencyGraph>::vertex_descriptor;
  using bimap_type       = boost::bimap<Vertex, unsigned>;

  bimap_type indices;
  std::map<Vertex, VertexDescriptor> vdm;

  unsigned vertexIndex = 0;

  for( auto it = S.begin_dimension(); it != S.end_dimension(); ++it )
  {
    // Vertices
    if( it->dimension() == 0 )
    {
      indices.insert( typename bimap_type::value_type( *it->begin(), vertexIndex++ ) );

      auto vertex             = boost::add_vertex( adjacencyGraph );
      vdm[ *it->begin() ]     = vertex;
    }

    // Edges
    else if( it->dimension() == 1 )
    {
      auto&& edge = *it;

      boost::add_edge( vdm.at( *( edge.begin() ) ),
                       vdm.at( *( edge.begin() + 1 ) ),
                       adjacencyGraph );
    }
  }

  return std::make_pair( indices, adjacencyGraph );
}

} // namespace detail


} // namespace aleph

#endif
