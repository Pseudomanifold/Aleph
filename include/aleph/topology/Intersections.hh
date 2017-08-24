#ifndef ALEPH_TOPOLOGY_INTERSECTIONS_HH__
#define ALEPH_TOPOLOGY_INTERSECTIONS_HH__

#include <aleph/math/Combinations.hh>

#include <algorithm>
#include <iterator>
#include <set>
#include <type_traits>
#include <vector>

namespace aleph
{

namespace topology
{

/**
  Intersects two simplices with each other. If the intersection is
  empty, the function returns the empty simplex. The data value of
  the simplex is left unspecified. It should be set by the client.
*/

template <class Simplex> Simplex intersect( const Simplex& s, const Simplex& t ) noexcept
{
  using VertexType = typename Simplex::VertexType;

  std::set<VertexType> sVertices( s.begin(), s.end() );
  std::set<VertexType> tVertices( t.begin(), t.end() );

  std::vector<VertexType> stVertices;
  stVertices.reserve( std::min( s.size(), t.size() ) );

  std::set_intersection( sVertices.begin(), sVertices.end(),
                         tVertices.begin(), tVertices.end(),
                         std::back_inserter( stVertices ) );

  return Simplex( stVertices.begin(), stVertices.end() );
}

/**
  Intersects a simplex with a simplicial complex. To this end, the
  intersection between all simplices in the simplicial complex and
  the current simplex will be calculated.

  The result of this intersection is a set of simplices.
*/

template <class SimplicialComplex, class Simplex> std::set<Simplex> intersect( const SimplicialComplex& K, const Simplex& s ) noexcept
{
  static_assert( std::is_same<Simplex, typename SimplicialComplex::ValueType>::value,
                 "Simplex type and simplicial complex value type must coincide" );

  // Shortcut: if the simplex is contained in the simplicial complex, we
  // report this as an intersection.
  if( K.contains(s) )
    return { *K.find(s) };

  std::set<Simplex> simplices;

  // Assuming that the simplicial complex is not malformed, it makes no
  // sense to check for intersections with simplices whose dimension is
  // larger than the input simplex.
  auto range = K.range( [  ] ( std::size_t /* d */ ) { return true; },
                        [&s] ( std::size_t d       ) { return d <= s.dimension(); } );

  for( auto it = range.first; it != range.second; ++it )
  {
    auto&& t = *it;
    auto u   = intersect(s,t);

    // Only insert non-empty intersections
    if( u )
      simplices.insert( u );
  }

  return simplices;
}

/**
  Intersects a simplex with a simplicial complex while constraining the
  dimensionality of the intersection. The idea is that a 0-simplex will
  only lead to 0-dimensional intersections, so there is no need to take
  a look at higher-dimensional simplices.

  Note that this assumes that the given simplicial complex is a complex
  in the sense of mathematics---it *must* contain every face of *every*
  simplex.

  The result of this intersection is a set of simplices.
*/

template <class SimplicialComplex, class Simplex> std::set<Simplex> intersectWithConstrainedDimension( const SimplicialComplex& K, const Simplex& s ) noexcept
{
  static_assert( std::is_same<Simplex, typename SimplicialComplex::ValueType>::value,
                 "Simplex type and simplicial complex value type must coincide" );

  std::set<Simplex> simplices;

  auto range = K.range( s.dimension() );

  for( auto it = range.first; it != range.second; ++it )
  {
    auto&& t = *it;
    auto u   = intersect(s,t);

    // Only insert non-empty intersections
    if( u )
      simplices.insert( u );
  }

  return simplices;
}

/**
  Calculates the *last* lexicographical intersection of the given
  simplex with the given simplicial complex. *Last* refers to the
  fact that the function attempts to find the intersection of the
  highest possible dimensionality.

  @param K Simplicial complex
  @param s Simplex

  @returns Intersection or an empty simplex if there is none
*/

template <class SimplicialComplex, class Simplex> Simplex lastLexicographicalIntersection( const SimplicialComplex& K, const Simplex& s ) noexcept
{
  using Vertex         = typename Simplex::VertexType;
  using Vector         = std::vector<Vertex>;

  using Iterator       = typename Vector::const_iterator;
  using DifferenceType = typename Vector::difference_type;

  Vector vertices( s.begin(), s.end() );
  Simplex result;

  // The basic idea of this function is that *all* possible subsets of
  // the simplex need to be evaluated, as only one of them is possibly
  // a match for the intersection.
  //
  // If the simplicial complex is large and the dimension of the input
  // simplex is small, it is much cheaper to query the complex for the
  // simplex rather than actually *calculating* intersections manually
  // with all simplices of the complex.
  for( std::size_t d = s.size(); d >= 1; d-- )
  {
    aleph::math::for_each_combination( vertices.begin(), vertices.begin() + DifferenceType( d ), vertices.end(),
      [&K, &result] ( Iterator first, Iterator last )
      {
        auto pos = K.find( Simplex( first, last ) );
        if( pos != K.end() )
        {
          result = *pos;
          return true;
        }

        return false;
      }
    );

    // We may break the loop as soon as a non-empty simplex has been
    // identified.
    if( result )
      return result;
  }

  return result;
}

} // namespace topology

} // namespace aleph

#endif
