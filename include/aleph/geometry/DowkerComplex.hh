#ifndef ALEPH_GEOMETRY_DOWKER_COMPLEX_HH__
#define ALEPH_GEOMETRY_DOWKER_COMPLEX_HH__

#include <aleph/math/Combinations.hh>

#include <aleph/topology/filtrations/Data.hh>

#include <aleph/topology/Simplex.hh>
#include <aleph/topology/SimplicialComplex.hh>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>

#include <boost/graph/floyd_warshall_shortest.hpp>
#include <boost/graph/johnson_all_pairs_shortest.hpp>

#include <algorithm>
#include <iterator>
#include <limits>
#include <unordered_map>
#include <vector>

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

template <class T, class I = std::size_t> struct Pair
{
  I p; // first index
  I q; // second index
  T w; // weight
};

template <class D, class V> struct Vertex
{
  V p; // vertex index
  D w; // weight
};

template <class D, class V> void swap( Vertex<D, V>& first, Vertex<D, V>& second )
{
  std::swap( first.p, second.p );
  std::swap( first.w, second.w );
}

} // namespace detail

/**
  Calculates a set of admissible pairs from a matrix of weights and
  a given distance threshold. The matrix of weights does *not* have
  to satisfy symmetry constraints.

  @param W Weighted adjacency matrix
  @param R Maximum weight
*/

template <class Matrix, class T> std::vector< detail::Pair<T> > admissiblePairs( const Matrix& W, T R )
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

  std::vector< Pair<T> > pairs;

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
        pairs.push_back( {i, j, static_cast<T>( D[i][j] ) } );
    }
  }

  return pairs;
}

/**
  Creates a Dowker sink complex and a Dowker source complex from a given
  set of admissible pairs. A *general* Dowker complex contains a simplex
  if all of its vertices satisfy the admissibility condition.

  @param pairs     Set of admissible pairs
  @param dimension Maximum dimension for expansion. If set to zero, will
                   expand the complex to its maximum dimension.
*/

template <class V, class D, class T>
std::pair<
  topology::SimplicialComplex< topology::Simplex<D, V> >,
  topology::SimplicialComplex< topology::Simplex<D, V> >
> buildDowkerSinkSourceComplexes( const std::vector<detail::Pair<T> >& pairs,
                                  unsigned dimension = 0 )
{
  using namespace detail;

  using Simplex           = topology::Simplex<D, V>;
  using SimplicialComplex = topology::SimplicialComplex<Simplex>;

  using VertexType     = V;
  VertexType maxVertex = VertexType();

  for( auto&& pair : pairs )
  {
    maxVertex = std::max(maxVertex, VertexType(pair.p) );
    maxVertex = std::max(maxVertex, VertexType(pair.q) );
  }

  using Vertex = Vertex<D, V>;

  // Keep track of the mapping induced by fixing either the source
  // points or the sink points.
  std::unordered_map< VertexType, std::vector<Vertex> > sourceBasePointMap;
  std::unordered_map< VertexType, std::vector<Vertex> > sinkBasePointMap;

  for( auto&& pair : pairs )
  {
    auto&& p = pair.p;
    auto&& q = pair.q;

    sourceBasePointMap[ VertexType(p) ].push_back( { VertexType(q), pair.w } );
    sinkBasePointMap[ VertexType(q) ].push_back( { VertexType(p), pair.w } );
  }

  // Auxiliary weight calculation lambda function ----------------------
  //
  // This function calculates the *maximum* weight of a range of
  // vertices. It is used to determine the weight of a simplex.

  using Iterator = typename std::vector<Vertex>::const_iterator;
  auto getWeight = [] ( Iterator begin, Iterator end )
  {
    using DataType  = D;
    DataType weight = std::numeric_limits<DataType>::lowest();

    for( auto it = begin; it != end; ++it )
      weight = std::max( weight, it->w );

    return weight;
  };

  auto makeSimplices = [&dimension, &getWeight] ( const std::unordered_map< VertexType, std::vector<Vertex> >& map )
  {
    std::vector<Simplex> simplices;
    std::unordered_map<Simplex, D> simplex_to_weight;

    for( auto&& pair : map )
    {
      auto vertices            = pair.second;
      std::size_t maxDimension = 0;

      if( dimension == 0 )
        maxDimension = vertices.size();
      else
        maxDimension = dimension + 1;

      for( std::size_t d = std::min( vertices.size(), maxDimension ); d >= 1; d-- )
      {
        math::for_each_combination( vertices.begin(), vertices.begin() + d, vertices.end(),
          [&simplex_to_weight, &getWeight] ( Iterator first, Iterator last )
          {
            std::vector<V> vertices_;
            vertices_.reserve( typename std::vector<V>::size_type( std::distance( first, last ) ) );

            for( auto it = first; it != last; ++it )
              vertices_.push_back( it->p );

            Simplex s( vertices_.begin(), vertices_.end() );

            if( simplex_to_weight.find(s) == simplex_to_weight.end() )
              simplex_to_weight[s] = getWeight(first, last);
            else
              simplex_to_weight[s] = std::min( simplex_to_weight[s], getWeight(first, last) );

            return false;
          }
        );
      }
    }

    for( auto&& pair : simplex_to_weight )
    {
      auto s = pair.first;
      s.setData( pair.second );

      simplices.push_back( s );
    }

    return simplices;
  };

  auto sourceEdges = makeSimplices( sourceBasePointMap );
  auto sinkEdges   = makeSimplices( sinkBasePointMap );

  SimplicialComplex dowkerSourceComplex( sourceEdges.begin(), sourceEdges.end() );
  SimplicialComplex dowkerSinkComplex  ( sinkEdges.begin()  , sinkEdges.end()   );

  dowkerSourceComplex.sort( topology::filtrations::Data<Simplex>() );
  dowkerSinkComplex.sort( topology::filtrations::Data<Simplex>() );

  return std::make_pair( dowkerSourceComplex, dowkerSinkComplex );
}

} // namespace geometry

} // namespace aleph

#endif
