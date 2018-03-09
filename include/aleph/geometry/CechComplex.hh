#ifndef ALEPH_GEOMETRY_CECH_COMPLEX_HH__
#define ALEPH_GEOMETRY_CECH_COMPLEX_HH__

#include <aleph/geometry/RipsExpander.hh>
#include <aleph/geometry/RipsSkeleton.hh>

#include <aleph/math/Combinations.hh>

#include <aleph/topology/Simplex.hh>
#include <aleph/topology/SimplicialComplex.hh>

#include <aleph/topology/filtrations/Data.hh>

#include <algorithm>
#include <vector>

namespace aleph
{

namespace geometry
{

template <class NearestNeighbours> auto buildCechComplex3D(
  const NearestNeighbours& nn,
  typename NearestNeighbours::ElementType epsilon,
  unsigned dimension ) -> topology::SimplicialComplex< topology::Simplex<typename NearestNeighbours::ElementType, typename NearestNeighbours::IndexType> >
{
  using ElementType       = typename NearestNeighbours::ElementType;
  using IndexType         = typename NearestNeighbours::IndexType;
  using Simplex           = topology::Simplex<ElementType, IndexType>;
  using SimplicialComplex = topology::SimplicialComplex<Simplex>;

  geometry::RipsSkeleton<NearestNeighbours> ripsSkeleton;

  auto skeleton
    = ripsSkeleton( nn, epsilon );

  // Set up vertices for a combinatorial search over *all* potential
  // simplices.
  std::vector<IndexType> vertices( nn.size() );
  std::iota( vertices.begin(), vertices.end(),
             IndexType(0) );

  using DifferenceType = typename decltype(vertices)::difference_type;
  using Iterator       = typename decltype(vertices)::const_iterator;

  std::vector<Simplex> simplices;

  math::for_each_combination( vertices.begin(), vertices.begin() + DifferenceType(3); vertices.end(),
    [&simplices] ( Iterator first, Iterator last )
    {
      Simplex s( first, last );

      simplices.push_back( s );
      return false;
    }
  );

  return SimplicialComplex( simplices.begin(), simplices.end() );
}

} // namespace geometry

} // namespace aleph


#endif
