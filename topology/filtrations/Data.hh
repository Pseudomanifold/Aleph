#ifndef ALEPH_TOPOLOGY_FILTRATIONS_DATA_HH__
#define ALEPH_TOPOLOGY_FILTRATIONS_DATA_HH__

#include <functional>

namespace aleph
{

namespace filtrations
{

template <
  class Simplex,
  class Compare = std::less<typename Simplex::DataType>
> class Data
{
public:
  bool operator()( const Simplex& s, const Simplex& t ) const
  {
    if( s.data() == t.data() )
    {
      // Default to lexicographical comparison if the two dimensions
      // coincide. We do not have any other choice here.
      if( s.dimension() == t.dimension() )
        return s < t;

      // Faces need to precede cofaces in order to obtain a valid
      // filtration.
      else
        return s.dimension() < t.dimension();
    }
    else
      return Compare()( s.data(), t.data() );
  }
};

}

}

#endif
