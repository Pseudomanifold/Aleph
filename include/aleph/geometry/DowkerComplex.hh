#ifndef ALEPH_GEOMETRY_DOWKER_COMPLEX_HH__
#define ALEPH_GEOMETRY_DOWKER_COMPLEX_HH__

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>

#include <vector>

// FIXME: remove after debugging
#include <iostream>
#include <boost/graph/graphviz.hpp>

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

  (void) R;

  boost::write_graphviz( std::cout, G );
  return {};
}

} // namespace geometry

} // namespace aleph

#endif
