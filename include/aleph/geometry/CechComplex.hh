#ifndef ALEPH_GEOMETRY_CECH_COMPLEX_HH__
#define ALEPH_GEOMETRY_CECH_COMPLEX_HH__

#include <aleph/math/Combinations.hh>

#include <aleph/topology/Simplex.hh>
#include <aleph/topology/SimplicialComplex.hh>

#include <aleph/topology/filtrations/Data.hh>

#include <aleph/external/Miniball.hpp>

#include <numeric>
#include <vector>

#include <cmath>

namespace aleph
{

namespace geometry
{

template <class Container> auto buildCechComplex( const Container& container, typename Container::ElementType r ) -> topology::SimplicialComplex< topology::Simplex<typename Container::ElementType, typename Container::IndexType> >
{
  using ElementType       = typename Container::ElementType;
  using IndexType         = typename Container::IndexType;
  using Simplex           = topology::Simplex<ElementType, IndexType>;
  using SimplicialComplex = topology::SimplicialComplex<Simplex>;

  // Set up vertices for a combinatorial search over *all* potential
  // simplices.
  std::vector<IndexType> vertices( container.size() );
  std::iota( vertices.begin(), vertices.end(),
             IndexType(0) );

  std::vector<Simplex> simplices;

  // Add 0-skeleton ----------------------------------------------------

  for( auto&& vertex : vertices )
    simplices.push_back( Simplex( vertex ) );

  // Add 1-skeleton ----------------------------------------------------
  //
  // We could obtain this faster by using nearest-neighbour queries, but
  // let's remain in the same domain and calculate the complex similarly
  // in *every* dimension.

  auto D               = container.dimension();
  auto n               = container.size();
  auto R               = r * r;
  using DifferenceType = typename decltype(vertices)::difference_type;
  using Iterator       = typename decltype(vertices)::const_iterator;

  for( unsigned d = 2; d <= n; d++ )
  {
    auto previous_n = simplices.size();

    math::for_each_combination( vertices.begin(), vertices.begin() + DifferenceType(d), vertices.end(),
      [&container, &simplices, &D, &R] ( Iterator first, Iterator last )
      {
        std::vector< std::vector<ElementType> > points;
        for( Iterator it = first; it != last; ++it )
          points.push_back( container[ *it ] );

        using PointIterator      = typename decltype(points)::const_iterator;
        using CoordinateIterator = typename std::vector<ElementType>::const_iterator;
        using Miniball           = Miniball::Miniball< Miniball::CoordAccessor<PointIterator, CoordinateIterator> >;

        Miniball mb( static_cast<int>(D), points.begin(), points.end() );
        if( mb.squared_radius() <= R )
        {
          Simplex s( first, last, ElementType( 2 * std::sqrt( mb.squared_radius() ) ) );
          simplices.push_back( s );
        }

        return false;
      }
    );

    // Check whether we can stop adding simplices. If no additional
    // simplices have been created during the last iteration, there
    // is no need to continue the iteration.
    auto current_n = simplices.size();
    if( previous_n == current_n )
      break;
  }

  SimplicialComplex K( simplices.begin(), simplices.end() );
  K.sort( topology::filtrations::Data<Simplex>() );

  return K;
}

} // namespace geometry

} // namespace aleph


#endif
