#ifndef ALEPH_PERSISTENCE_PAIRING_HH__
#define ALEPH_PERSISTENCE_PAIRING_HH__

#include <algorithm>
#include <limits>
#include <utility>
#include <vector>

namespace aleph
{

/**
  @class PersistencePairing
  @brief Container for index-based persistence pairings

  This class is a general-purpose container for pairings based on
  persistent homology. It consists of pairs of indices that refer
  to the paired simplices (or critical points) calculated using a
  persistent homology algorithm, for example.

  The class is purposefully kept simple and represents 'unpaired'
  simplices using a very large value. More precisely,

      std::numeric_limits<Index>::max()

  is used.
*/

template <class Index> class PersistencePairing
{
public:

  // Typedefs & aliases ------------------------------------------------

  using IndexType     = Index;

  using ValueType     = std::pair<Index, Index>;
  using ContainerType = std::vector<ValueType>;

  using ConstIterator = typename ContainerType::const_iterator;
  using Iterator      = typename ContainerType::iterator;

  using value_type    = ValueType;

  // Iterators ---------------------------------------------------------

  ConstIterator begin() const { return _pairs.begin(); }
  Iterator      begin()       { return _pairs.begin(); }

  ConstIterator end() const   { return _pairs.end(); }
  Iterator      end()         { return _pairs.end(); }

  // Container modification --------------------------------------------

  void add( Index birth )
  {
    _pairs.push_back( std::make_pair( birth, std::numeric_limits<Index>::max() ) );
  }

  void add( Index birth, Index destruction )
  {
    _pairs.push_back( std::make_pair( birth, destruction ) );
  }

  Iterator erase( Iterator position )
  {
    return _pairs.erase( position );
  }

  Iterator erase( Iterator begin, Iterator end )
  {
    return _pairs.erase( begin, end );
  }

  // Comparison operators ----------------------------------------------

  bool operator==( const PersistencePairing<Index>& other ) const
  {
    return _pairs == other._pairs;
  }

  bool operator!=( const PersistencePairing<Index>& other ) const
  {
    return !( this->operator==( other ) );
  }

  // Queries -----------------------------------------------------------

  /**
    Returns the iterator corresponding to a given pair of indices. If
    the pairing does not contain the pair the function returns an end
    iterator.
  */

  ConstIterator find( Index creator, Index destroyer ) const noexcept
  {
    return std::find( _pairs.begin(), _pairs.end(),
                      std::make_pair( creator, destroyer ) );
  }

  /**
    Returns the iterator corresponding to a given creator index. If the
    creator cannot be found in the pairing, the function returns an end
    iterator.
  */

  ConstIterator find( Index creator ) const noexcept
  {
    return std::find_if( _pairs.begin(), _pairs.end(),
                         [&creator] ( const ValueType& pair )
                         {
                           return pair.first == creator;
                         } );

  }

  bool contains( Index creator, Index destroyer ) const
  {
    return this->find( creator, destroyer ) != _pairs.end();
  }

  bool contains( Index creator ) const
  {
    return this->find( creator ) != _pairs.end();
  }

  std::size_t size() const
  {
    return _pairs.size();
  }

  bool empty() const
  {
    return _pairs.empty();
  }

private:
  ContainerType _pairs;
};

} // namespace aleph

#endif
