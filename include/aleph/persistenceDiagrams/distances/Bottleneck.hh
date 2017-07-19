#ifndef ALEPH_PERSISTENCE_DIAGRAMS_DISTANCES_BOTTLENECK_HH__
#define ALEPH_PERSISTENCE_DIAGRAMS_DISTANCES_BOTTLENECK_HH__

#include <aleph/geometry/distances/Infinity.hh>

#include <aleph/persistenceDiagrams/PersistenceDiagram.hh>
#include <aleph/persistenceDiagrams/distances/detail/Orthogonal.hh>

#include <boost/iterator/counting_iterator.hpp>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/max_cardinality_matching.hpp>

namespace aleph
{

namespace distances
{

namespace detail
{

template <class T> struct Edge
{
  std::size_t source;
  std::size_t target;

  T weight;

  Edge( std::size_t s, std::size_t t, T w )
    : source( s )
    , target( t )
    , weight( w )
  {
  }

  bool operator<( const Edge& other ) const
  {
    return weight < other.weight;
  }
};

template <class T> struct CheckMatchingCardinality
{
  using GraphType          = boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS>;
  using MatchingVectorType = std::vector<boost::graph_traits<GraphType>::vertex_descriptor>;

  CheckMatchingCardinality( std::size_t size, typename std::vector<Edge<T> >::const_iterator begin )
    : _maximumSize( size )
    , _last( begin )
    , _graph( 2 * size )
    , _mates( 2 * size )
  {
    boost::add_edge( begin->source, begin->target, _graph );
  }

  bool operator()( typename std::vector<Edge<T> >::const_iterator /* it1 */, typename std::vector<Edge<T> >::const_iterator it2 )
  {
    // The new edge lies beyond the edges that are already known, so the
    // edges between the last edge and the new iterator position need to
    // be added.
    if( it2 > _last )
    {
      do
      {
        ++_last;
        boost::add_edge( _last->source, _last->target, _graph );
      }
      while( _last != it2 );
    }

    // The new edge lies behind the edges that are already known, so the
    // surplus edges need to be removed.
    else
    {
      do
      {
        boost::remove_edge( _last->source, _last->target, _graph );
        --_last;
      }
      while( _last != it2 );
    }

    boost::edmonds_maximum_cardinality_matching( _graph,
                                                 &_mates[0] );

    // Look out for _perfect matchings_ in the bipartite graph. Any other
    // maximum cardinality matching does not qualify for the Bottleneck
    // distance.
    return boost::matching_size( _graph, &_mates[0] ) == _maximumSize;
  }

  std::size_t _maximumSize;

  typename std::vector<Edge<T> >::const_iterator _last;

  GraphType          _graph; // Input data
  MatchingVectorType _mates; // Edges of the matching
};

} // namespace detail

/**
  Calculates the Bottleneck distance between two persistence diagrams.
  The algorithm used for this involves checking a (complete) bipartite
  graph for perfect matchings.

  A brief description of the algoritmh is given in

    Computational Topology
    Herbert Edelsbrunner and John Harer

  on page 191.

  The implementation has been inspired by Dmitriy Morozov's "Dionysus"
  framework.

  @param D1 First persistence diagram
  @param D2 Second persistence diagram

  @returns Bottleneck distance between the two persistence diagrams
*/


template <
  class DataType,
  class Distance = InfinityDistance<DataType>
> DataType bottleneckDistance( const PersistenceDiagram<DataType>& D1,
                               const PersistenceDiagram<DataType>& D2 )
{
  auto n           = D1.size();
  auto m           = D2.size();
  auto maximumSize = n + m;

  using SizeType = decltype(n);
  using Edge     = detail::Edge<DataType>;

  std::vector<Edge> edges;

  // Diagonal edges ----------------------------------------------------

  for( SizeType i = n; i < maximumSize; i++ )
    for( SizeType j = maximumSize + i; j < 2 * maximumSize; j++ )
      edges.push_back( Edge( static_cast<std::size_t>(i), static_cast<std::size_t>(j), DataType() ) );

  SizeType i = 0;

  // Edges between regular points --------------------------------------

  Distance distance;

  for( auto it1 = D1.begin(); it1 != D1.end(); ++it1 )
  {
    auto j = maximumSize;

    for( auto it2 = D2.begin(); it2 != D2.end(); ++it2 )
    {
      auto weight = distance( *it1, *it2 );

      edges.push_back( Edge( static_cast<std::size_t>(i), static_cast<std::size_t>(j),
                             weight ) );

      ++j;
    }

    ++i;
  }

  // Edges between points and their projections ------------------------

  i = 0;

  for( auto it1 = D1.begin(); it1 != D1.end(); ++it1 )
  {
    edges.push_back( Edge( static_cast<std::size_t>(i), static_cast<std::size_t>(maximumSize + m + i),
                           aleph::distances::detail::orthogonalDistance<Distance>( *it1 ) ) );

    ++i;
  }

  i = maximumSize;

  for( auto it2 = D2.begin(); it2 != D2.end(); ++it2 )
  {
    edges.push_back( Edge( static_cast<std::size_t>(n + i - maximumSize), static_cast<std::size_t>(i),
                           aleph::distances::detail::orthogonalDistance<Distance>( *it2 ) ) );
  }

  // Identify matchings ------------------------------------------------

  std::sort( edges.begin(), edges.end() );

  // Perform binary search over edge sets. Starting from the empty graph, use
  // more and more edges to find out the first graph that permits a maximum
  // cardinality matching.

  using EdgeIteratorType     = typename std::vector<Edge>::const_iterator;
  using CountingIteratorType = boost::counting_iterator<EdgeIteratorType>;

  auto itEdge = std::upper_bound( CountingIteratorType( edges.begin() ),
                                  CountingIteratorType( edges.end() ),
                                  edges.begin(),
                                  detail::CheckMatchingCardinality<DataType>( maximumSize, edges.begin() ) );

  return (*itEdge)->weight;
}

} // namespace distances

} // namespace aleph

#endif
