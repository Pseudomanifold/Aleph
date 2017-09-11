#ifndef ALEPH_TOPOLOGY_COMBINATORIAL_CURVATURE_HH__
#define ALEPH_TOPOLOGY_COMBINATORIAL_CURVATURE_HH__

#include <aleph/math/KahanSummation.hh>

#include <algorithm>
#include <iterator>
#include <set>
#include <vector>

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

template <class SimplicialComplex> auto commonCofaces(
  const SimplicialComplex& K,
  const typename SimplicialComplex::ValueType& s,
  const typename SimplicialComplex::ValueType& t )
  -> std::vector<typename SimplicialComplex::ValueType>
{
  if( s.dimension() != t.dimension() )
    return {};

  using Simplex = typename SimplicialComplex::ValueType;

  std::vector<Simplex> cofaces;

  auto range = K.range( s.dimension() + 1 );
  for( auto&& it = range.first; it != range.second; ++it )
  {
    if(    std::find( it->begin_boundary(), it->end_boundary(), s ) != it->end_boundary()
        && std::find( it->begin_boundary(), it->end_boundary(), t ) != it->end_boundary() )
    {
      cofaces.push_back( *it );
    }
  }

  return cofaces;
}

template <class SimplicialComplex> auto commonFaces(
  const SimplicialComplex& K,
  const typename SimplicialComplex::ValueType& s,
  const typename SimplicialComplex::ValueType& t )
  -> std::vector<typename SimplicialComplex::ValueType>
{
  if( s.dimension() != t.dimension() )
    return {};

  using Simplex = typename SimplicialComplex::ValueType;

  std::set<Simplex> commonBoundary;

  std::set_intersection( s.begin_boundary(), s.end_boundary(),
                         t.begin_boundary(), t.end_boundary(),
                         std::inserter( commonBoundary, commonBoundary.begin() )
  );

  std::vector<Simplex> commonFaces;

  for( auto&& simplex : commonBoundary )
    commonFaces.push_back( *K.find( simplex ) );

  return commonFaces;
}

template <class SimplicialComplex, class OutputIterator> void weightedCurvature( const SimplicialComplex& K,
                                                                                 OutputIterator result,
                                                                                 unsigned p = 1 )
{
  using Simplex    = typename SimplicialComplex::ValueType;
  using DataType   = typename Simplex::DataType;

  auto pair = K.range( p );

  for( auto&& itSimplex = pair.first; itSimplex != pair.second; ++itSimplex )
  {
    auto&& s = *itSimplex;

    // 1. Summand: Co-faces --------------------------------------------

    std::vector<DataType> weights_Cofaces;

    auto range = K.range( s.dimension() + 1 );
    for( auto&& itCoface = range.first; itCoface != range.second; ++itCoface )
    {
      if( std::find( itCoface->begin_boundary(), itCoface->end_boundary(), s ) != itCoface->end_boundary() )
        weights_Cofaces.push_back( s.data() / itCoface->data() );
    }

    // 2. Summand: Faces -----------------------------------------------

    std::vector<DataType> weights_Faces;

    for( auto itBoundary = s.begin_boundary(); itBoundary != s.end_boundary(); ++itBoundary )
    {
      auto itPosition = K.find( *itBoundary );
      if( itPosition != K.end() )
        weights_Faces.push_back( itPosition->data() / s.data() );
    }

    // 3. Summand: Parallel neighbours ---------------------------------

    range = K.range( s.dimension() );

    std::vector<DataType> weights_commonCofaces;
    std::vector<DataType> weights_commonFaces;

    for( auto itNeighbour = range.first; itNeighbour != range.second; ++itNeighbour )
    {
      if( *itNeighbour == s )
        continue;

      auto cofaces = commonCofaces( K, *itNeighbour, s );
      auto faces   = commonFaces  ( K, *itNeighbour, s );
      auto weight  = std::sqrt( s.data() * itNeighbour->data() );

      for( auto&& coface : cofaces )
        weights_commonCofaces.push_back( weight / coface.data() );

      for( auto&& face : faces )
        weights_commonFaces.push_back( face.data() / weight );
    }

    auto s11 = aleph::math::accumulate_kahan_sorted( weights_Cofaces.begin()      , weights_Cofaces.end()      , DataType() );
    auto s12 = aleph::math::accumulate_kahan_sorted( weights_Faces.begin()        , weights_Faces.end()        , DataType() );
    auto s21 = aleph::math::accumulate_kahan_sorted( weights_commonCofaces.begin(), weights_commonCofaces.end(), DataType() );
    auto s22 = aleph::math::accumulate_kahan_sorted( weights_commonFaces.begin()  , weights_commonFaces.end()  , DataType() );

    auto curvature = s.data() * ( (s11 + s12) - std::abs( s21 - s22 ) );
    *result++      = curvature;
  }
}

} // namespace topology

} // namespace aleph

#endif
