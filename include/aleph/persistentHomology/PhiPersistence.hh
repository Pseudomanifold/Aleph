#ifndef ALEPH_PERSISTENT_HOMOLOGY_PHI_PERSISTENCE_HH__
#define ALEPH_PERSISTENT_HOMOLOGY_PHI_PERSISTENCE_HH__

#include <aleph/topology/SimplicialComplex.hh>

#include <initializer_list>
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
      else if( _values[k] >= static_cast<int>( k-1 ) )
        _values[k] = static_cast<int>( k-1 );
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

} // namespace aleph

#endif
