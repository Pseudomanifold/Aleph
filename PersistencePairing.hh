#ifndef ALEPH_PERSISTENCE_PAIRING_HH__
#define ALEPH_PERSISTENCE_PAIRING_HH__

#include <limits>
#include <utility>
#include <vector>

namespace aleph
{

template <class Index> class PersistencePairing
{
public:

  // Typedefs & aliases ------------------------------------------------

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

  // Queries -----------------------------------------------------------

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
