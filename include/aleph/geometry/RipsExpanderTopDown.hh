#ifndef ALEPH_GEOMETRY_RIPS_EXPANDER_TOP_DOWN_HH__
#define ALEPH_GEOMETRY_RIPS_EXPANDER_TOP_DOWN_HH__

#include <algorithm>
#include <limits>
#include <list>
#include <vector>

#include "topology/MaximalCliques.hh"

namespace aleph
{

namespace geometry
{

namespace detail
{

/*
  Originally written by Thomas Draper for an article of Mark Nelson in the
  C/C++ users journal. See

    http://marknelson.us/2002/03/01/next-permutation/

  for more information.
*/

template <class Iterator> bool next_combination( const Iterator first, Iterator k, const Iterator last )
{
  if( first == last || first == k || last == k )
    return false;

  auto itr1 = first;
  auto itr2 = last;
  ++itr1;

  if(last == itr1)
    return false;

  itr1 = last;
  --itr1;
  itr1 = k;
  --itr2;

  while( first != itr1 )
  {
    if( *--itr1 < *itr2 )
    {
      auto j = k;
      while( !(*itr1 < *j) )
        ++j;

      std::iter_swap( itr1, j );

      ++itr1;
      ++j;
      itr2 = k;

      std::rotate( itr1, j, last );

      while( last != j )
      {
        ++j;
        ++itr2;
      }

      std::rotate( k, itr2, last );
      return true;
    }
  }

  std::rotate( first, k, last );
  return false;
}

} // namespace detail

template <class SimplicialComplex> class RipsExpanderTopDown
{
public:
  using Simplex           = typename SimplicialComplex::ValueType;
  using DataType          = typename Simplex::DataType;
  using VertexType        = typename Simplex::VertexType;

  SimplicialComplex operator()( const SimplicialComplex& K, unsigned kMax )
  {
    return this->operator()( K, kMax, 0 );
  }

  SimplicialComplex operator()( const SimplicialComplex& K, unsigned kMax, unsigned kMin )
  {
    auto maximalCliques = aleph::topology::maximalCliquesKoch( K );

    std::list<Simplex> simplices;

    for( auto&& clique : maximalCliques )
    {
      auto C = std::vector<VertexType>( clique.begin(), clique.end() );

      for( unsigned k = kMin + 1; k <= std::min( kMax + 1, unsigned( C.size() ) ); k++ )
      {
        do
        {
          simplices.push_back( Simplex( C.begin(), C.begin() + k ) );
        }
        while( detail::next_combination( C.begin(), C.begin() + k, C.end() ) );
      }
    }

    SimplicialComplex S;
    S.insert( simplices.begin(), simplices.end() );
    return S;
  }

  // Weight assignment -------------------------------------------------

  /**
    Given a simplicial complex K, uses another simplicial complex S to
    look up the weights that are to be assigned to K. The look up uses
    simplices from the specified dimension only and uses their maximum
    weight as the weight for the simplex in K.
  */

  SimplicialComplex assignMaximumWeight( const SimplicialComplex& K, const SimplicialComplex& S, unsigned dimension = 1 )
  {
    SimplicialComplex R;

    for( auto s : K )
    {
      std::vector<VertexType> vertices( s.begin(), s.end() );
      std::sort( vertices.begin(), vertices.end() );

      DataType weight = std::numeric_limits<DataType>::lowest();

      // Unable to assign a 0-simplex a weight using nothing but the
      // 1-simplices, for example.
      if( s.dimension() < dimension )
      {
        R.push_back( s );
        continue;
      }

      do
      {
        auto lookup = Simplex( vertices.begin(), vertices.begin() + dimension + 1 );
        auto it     = S.find( lookup );

        if( it != S.end() )
          weight = std::max( weight, it->data() );
      }
      while( detail::next_combination( vertices.begin(), vertices.begin() + dimension + 1, vertices.end() ) );

      s.setData( weight );
      R.push_back( s );
    }

    return R;
  }
};

} // namespace geometry

} // namespace aleph

#endif
