#ifndef ALEPH_TOPOLOGY_COMBINATORIAL_CURVATURE_HH__
#define ALEPH_TOPOLOGY_COMBINATORIAL_CURVATURE_HH__

#include <algorithm>
#include <iterator>
#include <set>

namespace aleph
{

namespace topology
{

template <class Simplex> bool hasFace( const Simplex& s, const Simplex& f )
{
  return std::find( s.begin_boundary(), s.end_boundary(), f ) != s.end_boundary();
}

template <class SimplicialComplex> bool parallelNeighbours( const SimplicialComplex& K,
                                                            const typename SimplicialComplex::ValueType& s,
                                                            const typename SimplicialComplex::ValueType& t )
{
  if( s.dimension() != t.dimension() )
    return false;

  using Simplex = typename SimplicialComplex::ValueType;

  std::set<Simplex> commonBoundary;

  std::set_intersection( s.begin_boundary(), s.end_boundary(),
                         t.begin_boundary(), t.end_boundary(),
                         std::inserter( commonBoundary, commonBoundary.begin() )
  );

  bool condition1 = !commonBoundary.empty();
  bool condition2 = false;

  auto iterators = K.range( s.dimension() + 1 );

  for( auto&& it = iterators.first; it != iterators.second; ++it )
  {
    if( hasFace( *it, s ) && hasFace( *it, t ) )
    {
      condition2 = true;
      break;
    }
  }

  // Only one of the conditions is allowed to be true for the two
  // simplices to be considered parallel neighbours.
  return condition1 != condition2;
}

template <class SimplicialComplex, class OutputIterator> void curvature( const SimplicialComplex& K,
                                                                         OutputIterator result,
                                                                         unsigned p = 1 )
{
  using Simplex    = typename SimplicialComplex::ValueType;
  using VertexType = typename Simplex::VertexType;
  auto pair        = K.range( p );

  for( auto&& it = pair.first; it != pair.second; ++it )
  {
    auto&& s                   = *it;
    auto isCoface              = [&s]    ( const Simplex& t ) { return hasFace( t, s ); };
    auto isParallelNeighbour   = [&K,&s] ( const Simplex& t ) { return parallelNeighbours( K, s, t ); };
    auto range                 = K.range( s.dimension() + 1 );
    auto numCofaces            = std::count_if( range.first, range.second, isCoface );
    range                      = K.range( s.dimension() );
    auto numParallelNeighbours = std::count_if( range.first, range.second, isParallelNeighbour );

    *result++ = VertexType(   numCofaces
                            + s.size()
                            - numParallelNeighbours );
  }
}

} // namespace topology

} // namespace aleph

#endif
