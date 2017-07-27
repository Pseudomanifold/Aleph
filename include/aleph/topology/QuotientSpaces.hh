#ifndef ALEPH_TOPOLOGY_QUOTIENT_SPACES_HH__
#define ALEPH_TOPOLOGY_QUOTIENT_SPACES_HH__

#include <iterator>
#include <set>
#include <vector>

namespace aleph
{

namespace topology
{

/**
  Calculates the cone over a simplicial complex. The function creates
  a new cone vertex that is guaranteed not to occur in the simplicial
  complex. An optional functor may be used to set the weight for each
  new simplex. Given a simplex $\sigma$, the functor should calculate
  the weight of $\sigma \cup \{v_c\}$, where $v_c$ is the cone vertex
  of the cone complex. The functor thus needs to return the weight of
  the cone vertex if presented with an empty simplex.

  @param K Simplicial complex
  @param f Functor for assigning weights (see above)

  @returns The cone over the given simplicial complex. Note that this
           operation may require re-sorting the complex.
*/

template <class SimplicialComplex, class Functor> SimplicialComplex cone( const SimplicialComplex& K, Functor f )
{
  using Simplex    = typename SimplicialComplex::ValueType;
  using VertexType = typename Simplex::VertexType;

  if( K.empty() )
    return {};

  VertexType coneVertex = VertexType();

  {
    std::set<VertexType> vertices;
    K.vertices( std::inserter( vertices, vertices.begin() ) );

    if( vertices.empty() )
      return {};

    coneVertex = static_cast<VertexType>( *vertices.rbegin() + 1 );
  }

  auto L = K;

  std::vector<Simplex> simplices;
  simplices.reserve( K.size() + 1 );

  simplices.emplace_back( Simplex( coneVertex, f( Simplex() ) ) );

  for( auto&& simplex : K )
  {
    std::vector<VertexType> vertices( simplex.begin(), simplex.end() );
    vertices.reserve( simplex.size() + 1 );
    vertices.emplace_back( coneVertex );

    simplices.emplace_back(
      Simplex( vertices.begin(), vertices.end(),
               f( simplex ) ) );
  }

  L.insert( simplices.begin(), simplices.end() );
  return L;
}

/**
  Overloaded variant for calculating the cone of a simplicial complex.
  This function assigns all new simplices of the cone a data value of
  zero.

  @see aleph::topology::cone()
*/

template <class SimplicialComplex> SimplicialComplex cone( const SimplicialComplex& K )
{
  using Simplex  = typename SimplicialComplex::ValueType;
  using DataType = typename Simplex::DataType;

  auto f = [] ( const Simplex& /* s */ )
              {
                return DataType();
              };

  return cone(K, f);
}

} // namespace topology

} // namespace alpeh

#endif
