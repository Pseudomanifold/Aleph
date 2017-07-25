#ifndef ALEPH_PERSISTENT_HOMOLOGY_PHI_PERSISTENCE_HH__
#define ALEPH_PERSISTENT_HOMOLOGY_PHI_PERSISTENCE_HH__

#include <aleph/persistenceDiagrams/PersistenceDiagram.hh>

#include <aleph/topology/Conversions.hh>
#include <aleph/topology/Intersections.hh>
#include <aleph/topology/SimplicialComplex.hh>

#include <aleph/persistentHomology/algorithms/Twist.hh>

#include <initializer_list>
#include <map>
#include <stdexcept>
#include <utility>
#include <vector>

namespace aleph
{

/**
  Partitions a simplicial complex according to its $\phi$-persistence
  values. This follows the persistent intersection homology algorithm
  in:

    Persistent Intersection Homology
    Paul Bendich and John Harer

  The function expects a simplicial complex $K$ and a function $\phi$
  that determines whether a simplex is proper or not. The function is
  going to create a new simplicial complex. This complex contains all
  proper simplices (in their original order) followed by all improper
  ones.
*/

template <class Simplex, class Function> std::pair<topology::SimplicialComplex<Simplex>, std::size_t> partition( const topology::SimplicialComplex<Simplex>& K, Function phi )
{
  topology::SimplicialComplex<Simplex> L;

  for( auto&& simplex : K )
  {
    if( phi(simplex) )
      L.push_back( simplex );
  }

  auto s = L.size();

  for( auto&& simplex : K )
  {
    if( !phi(simplex) )
      L.push_back( simplex );
  }

  return std::make_pair( L, s );
}

/**
  Models a perversity in the sense of intersection homology. The class
  ensures that all values satisfy

  \f
    -1 \leq p_k \leq k-1
  \f
*/

class Perversity
{
public:
  Perversity( std::initializer_list<int> values )
  {
    _values.assign( values.begin(), values.end() );

    for( std::size_t k = 0; k < _values.size(); k++ )
    {
      if( _values[k] < -1 )
        _values[k] = -1;
      // There is an index shift going on here: since $k$ runs from $0$
      // to $d-1$, there is no need to shift the upper bound.
      else if( _values[k] > static_cast<int>( k ) )
        _values[k] = static_cast<int>( k );
    }
  }

  /**
    Queries the perversity value in a given dimension $d$. Invalid
    dimension values only cause the function to return a zero.
  */

  int operator()( std::size_t d ) const noexcept
  {
    if( d < _values.size() )
      return _values[d];
    else
      return 0;
  }

private:
  std::vector<int> _values;
};

template <class Simplex> auto calculateIntersectionHomology( const aleph::topology::SimplicialComplex<Simplex>& K,
                                                             const std::vector< aleph::topology::SimplicialComplex<Simplex> >& X,
                                                             const Perversity& p ) -> std::vector< PersistenceDiagram<typename Simplex::DataType> >
{
  // 0. Check consistency of strata
  // 1. Create allowability function based on the dimensionality of the
  //    intersection of simplices with individual strata.
  // 2. Calculate $phi$-persistence
  // 3. Convert the result into a persistence diagram.

  // Check consistency of strata ---------------------------------------
  //
  // The maximum dimension of the strata has to match the dimension of
  // the simplicial complex.

  {
    std::size_t minDimension = K.dimension();
    std::size_t maxDimension = 0;

    for( auto&& stratum : X )
    {
      minDimension = std::min( minDimension, stratum.dimension() );
      maxDimension = std::max( maxDimension, stratum.dimension() );
    }

    if( maxDimension != K.dimension() )
      throw std::runtime_error( "Invalid stratification" );
  }

  // Check whether simplex is allowable --------------------------------

  std::map<Simplex, bool> phi;

  {
    auto d = K.dimension();

    for( auto&& s : K )
    {
      bool admissible = true;

      for( std::size_t k = 0; k < d; k++ )
      {
        // The notation follows Bendich and Harer, so $i$ is actually
        // referring to a dimension instead of an index. Beware!
        auto i            = s.dimension();
        auto intersection = aleph::topology::intersect( X.at( d - 1 - k ), s );
        auto dimension    = intersection.empty() ? -1 : static_cast<long>( intersection.rbegin()->dimension() );
        admissible        = admissible && static_cast<long>( dimension ) <= static_cast<long>( i - k + p(k) );
      }

      phi[s] = admissible;
    }
  }

  // Partition according to allowable simplices ------------------------

  aleph::topology::SimplicialComplex<Simplex> L;
  std::size_t s = 0;

  std::tie( L, s ) =
    aleph::partition( K, [&phi] ( const Simplex& s )
                         {
                           return phi.at(s);
                         } );

  // Calculate persistent intersection homology ------------------------

  auto boundaryMatrix = aleph::topology::makeBoundaryMatrix( L, s );
  using IndexType     = typename decltype(boundaryMatrix)::Index;

  aleph::persistentHomology::algorithms::Twist algorithm;
  algorithm( boundaryMatrix );

  using DataType = typename Simplex::DataType;

  // FIXME: support more persistence diagrams...
  aleph::PersistenceDiagram<DataType> diagram;

  for( IndexType i = 0; i < boundaryMatrix.getNumColumns(); i++ )
  {
    auto lowestOne = boundaryMatrix.getMaximumIndex( i );
    if( lowestOne.second && lowestOne.first <= s )
    {
      auto&& simplex = L[i];
      auto&& partner = L[lowestOne.first];

      diagram.add( simplex.data(), partner.data() );
    }
    else if( !lowestOne.second )
    {
      auto&& simplex = L[i];
      diagram.add( simplex.data() );
    }
  }

  return { diagram };
}

} // namespace aleph

#endif
