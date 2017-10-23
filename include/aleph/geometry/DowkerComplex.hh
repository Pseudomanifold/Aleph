#ifndef ALEPH_GEOMETRY_DOWKER_COMPLEX_HH__
#define ALEPH_GEOMETRY_DOWKER_COMPLEX_HH__

#include <aleph/topology/Simplex.hh>
#include <aleph/topology/SimplicialComplex.hh>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>

#include <boost/graph/floyd_warshall_shortest.hpp>
#include <boost/graph/johnson_all_pairs_shortest.hpp>

#include <unordered_map>
#include <vector>

// FIXME: remove after debugging
#include <iostream>

namespace aleph
{

namespace geometry
{

namespace detail
{

// FIXME: the weight type of an edge should be configurable as an
// additional template parameter
using EdgeWeightProperty = boost::property<boost::edge_weight_t, double>;

// FIXME: the data type of the graph should be configurable as an
// additional template parameter.
using Graph = boost::adjacency_list<
  boost::vecS,
  boost::vecS,
  boost::directedS,
  boost::no_property,
  EdgeWeightProperty
>;

using VertexDescriptor = boost::graph_traits<Graph>::vertex_descriptor;

// FIXME: the index type should be configurable as an additional
// template parameter.
using Pair = std::pair<std::size_t, std::size_t>;

} // namespace detail

/**
  Calculates a set of admissible pairs from a matrix of weights and
  a given distance threshold. The matrix of weights does *not* have
  to satisfy symmetry constraints.

  @param W Weighted adjacency matrix
  @param R Maximum weight
*/

template <class Matrix, class T> std::vector<detail::Pair> admissiblePairs( const Matrix& W, T R )
{
  using namespace detail;

  // Convert matrix into a graph ---------------------------------------

  auto n          = W.size();
  using IndexType = decltype(n);

  detail::Graph G( n );

  for( IndexType i = 0; i < n; i++ )
  {
    for( IndexType j = 0; j < n; j++ )
    {
      if( W[i][j] > 0 )
      {
        EdgeWeightProperty weight = W[i][j];

        boost::add_edge( VertexDescriptor(i), VertexDescriptor(j),
                         weight,
                         G );
      }
    }
  }

  double density
    = static_cast<double>( boost::num_edges(G) ) / static_cast<double>( boost::num_vertices(G) * ( boost::num_vertices(G) - 1 ) );

  // This 'pseudo-matrix' contains the completion of the weight function
  // specified by the input matrix.
  std::vector< std::vector<double> > D( boost::num_vertices(G),
                                        std::vector<double>( boost::num_vertices(G) ) );

  if( density >= 0.5 )
    boost::floyd_warshall_all_pairs_shortest_paths( G, D );
  else
    boost::johnson_all_pairs_shortest_paths( G, D );

  std::vector<Pair> pairs;

  // Create admissible pairs -------------------------------------------
  //
  // A pair is admissible if it satisfies a reachability property,
  // meaning that the induced graph distance permits to reach both
  // vertices under the specified distance threshold.

  for( IndexType i = 0; i < n; i++ )
  {
    for( IndexType j = 0; j < n; j++ )
    {
      if( D[i][j] <= R )
        pairs.push_back( std::make_pair(i,j) );
    }
  }

  return pairs;
}

/**
  Creates a Dowker sink complex and a Dowker source complex from a given
  set of admissible pairs. A *general* Dowker complex contains a simplex
  if all of its vertices satisfy the admissibility condition.

  At present, this function only calculates the one-dimensional skeleton
  of the complexes.

  @param pairs Set of admissible pairs
*/

template <class V, class D>
std::pair<
  topology::SimplicialComplex< topology::Simplex<D, V> >,
  topology::SimplicialComplex< topology::Simplex<D, V> >
> buildDowkerSinkSourceComplexes( const std::vector<detail::Pair>& pairs )
{
  using Simplex           = topology::Simplex<D, V>;
  using SimplicialComplex = topology::SimplicialComplex<Simplex>;

  using VertexType     = V;
  VertexType maxVertex = VertexType();

  for( auto&& pair : pairs )
  {
    maxVertex = std::max(maxVertex, pair.first );
    maxVertex = std::max(maxVertex, pair.second);
  }

  // Keep track of the mapping induces by fixing either the source
  // points or the sink points.
  std::unordered_map< VertexType, std::vector<VertexType> > sourceBasePointMap;
  std::unordered_map< VertexType, std::vector<VertexType> > sinkBasePointMap;

  for( auto&& pair : pairs )
  {
    auto&& p = pair.first;
    auto&& q = pair.second;

    sourceBasePointMap[ VertexType(p) ].push_back( VertexType(q) );
    sinkBasePointMap[ VertexType(q) ].push_back( VertexType(p) );
  }


  auto makeEdges = [] ( const std::unordered_map< VertexType, std::vector<VertexType> >& map )
  {
    std::vector<Simplex> edges;

    for( auto&& pair : map )
    {
      auto&& p        = pair.first;
      auto&& vertices = pair.second;

      // TODO: what about the weights?
      for( auto it1 = vertices.begin(); it1 = vertices.end(); ++it1 )
        for( auto it2 = it1 + 1; it2 = vertices.end(); ++it2 )
          edges.push_back( Simplex( { *it1, *it2 } ) );
    }

    return edges;
  };


  auto sourceEdges = makeEdges( sourceBasePointMap );
  auto sinkEdges   = makeEdges( sinkBasePointMap );

  SimplicialComplex dowkerSourceComplex( sourceEdges.begin(), sourceEdges.end() );
  SimplicialComplex dowkerSinkComplex  ( sinkEdges.begin()  , sinkEdges.end()   );

  return std::make_pair( dowkerSourceComplex, dowkerSinkComplex );
}

} // namespace geometry

} // namespace aleph

#endif
