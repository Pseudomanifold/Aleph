#ifndef ALEPH_PERSISTENCE_PAIRING_HH__
#define ALEPH_PERSISTENCE_PAIRING_HH__

#include <algorithm>
#include <limits>
#include <utility>
#include <vector>

namespace aleph
{

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

  bool contains( Index creator, Index destroyer ) const
  {
    return std::find( _pairs.begin(), _pairs.end(),
                      std::make_pair( creator, destroyer ) ) != _pairs.end();
  }

  bool contains( Index creator ) const
  {
    return std::find_if( _pairs.begin(), _pairs.end(),
                         [this, &creator] ( const ValueType& pair )
                         {
                           return pair.first == creator;
                         } ) != _pairs.end();
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

}

#endif
