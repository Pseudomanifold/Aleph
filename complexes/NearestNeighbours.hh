#ifndef ALEPH_COMPLEXES_NEAREST_NEIGHBOURS_HH__
#define ALEPH_COMPLEXES_NEAREST_NEIGHBOURS_HH__

#include <vector>

namespace aleph
{

namespace complexes
{

template <class Wrapper, class ElementType, class IndexType> class NearestNeighbours
{
public:
  void radiusSearch( ElementType radius,
                     std::vector< std::vector<IndexType> >& indices,
                     std::vector< std::vector<ElementType> >& distances )
  {
    static_cast<Wrapper&>( *this ).radiusSearch( radius,
                                                 indices,
                                                 distances );
  }

  std::size_t size() const noexcept
  {
    return static_cast<const Wrapper&>( *this ).size(); 
  }
};

}

}

#endif
