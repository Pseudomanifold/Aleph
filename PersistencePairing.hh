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

  void add( Index birth )
  {
    _pairs.push_back( std::make_pair( birth, std::numeric_limits<Index>::max() ) );
  }

  void add( Index birth, Index destruction )
  {
    _pairs.push_back( std::make_pair( birth, destruction ) );
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
  std::vector< std::pair<Index, Index> > _pairs;
};

}

#endif
