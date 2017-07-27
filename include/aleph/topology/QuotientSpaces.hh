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
  new simplex.

  Given a simplex \f$\sigma\f$, the functor should calculate the weight
  of \f$\sigma \cup \{v_c\}\f$, where \f$v_c\f$ denotes the cone vertex
  of the cone complex. As a consequence, the functor is needs to return
  the weight of the cone vertex if presented with an empty simplex.

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

/**
  Calculates the suspension of a simplicial complex. This will result
  in a sort of double cone over the simplicial complex. Just like the
  cone calculation, this function also uses a functor to get the data
  values of each new simplex.

  @param K Simplicial complex
  @param f Functor for assigning weights (see above)

  @returns The suspension over the given simplicial complex. Note that
           this operation may require re-sorting the complex.

  @see aleph::topology::cone()
*/

template <class SimplicialComplex, class Functor> SimplicialComplex suspension( const SimplicialComplex& K, Functor f )
{
  using Simplex    = typename SimplicialComplex::ValueType;
  using VertexType = typename Simplex::VertexType;

  if( K.empty() )
    return {};

  VertexType upperConeVertex = VertexType();
  VertexType lowerConeVertex = VertexType();

  {
    std::set<VertexType> vertices;
    K.vertices( std::inserter( vertices, vertices.begin() ) );

    if( vertices.empty() )
      return {};

    upperConeVertex = static_cast<VertexType>( *vertices.rbegin() + 1 );
    lowerConeVertex = static_cast<VertexType>( upperConeVertex + 1 );
  }

  auto L = K;

  std::vector<Simplex> simplices;
  simplices.reserve( 2 * K.size() + 2 );

  simplices.emplace_back( Simplex( upperConeVertex, f( Simplex() ) ) );
  simplices.emplace_back( Simplex( lowerConeVertex, f( Simplex() ) ) );

  for( auto&& simplex : K )
  {
    std::vector<VertexType> vertices( simplex.begin(), simplex.end() );
    vertices.reserve( simplex.size() + 1 );
    vertices.emplace_back( upperConeVertex );

    simplices.emplace_back(
      Simplex( vertices.begin(), vertices.end(),
               f( simplex ) ) );

    vertices.back() = lowerConeVertex;

    simplices.emplace_back(
      Simplex( vertices.begin(), vertices.end(),
               f( simplex ) ) );
  }

  L.insert( simplices.begin(), simplices.end() );
  return L;
}

/**
  Overloaded variant for calculating the suspension, the double cone, of
  a simplicial complex. This function assigns all new simplices a weight
  of zero.

  @see aleph::topology::suspension()
*/

template <class SimplicialComplex> SimplicialComplex suspension( const SimplicialComplex& K )
{
  using Simplex  = typename SimplicialComplex::ValueType;
  using DataType = typename Simplex::DataType;

  auto f = [] ( const Simplex& /* s */ )
              {
                return DataType();
              };

  return suspension(K, f);
}


} // namespace topology

} // namespace alpeh

#endif
