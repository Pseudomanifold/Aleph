#ifndef ALEPH_COMPLEXES_FLANN_HH__
#define ALEPH_COMPLEXES_FLANN_HH__

#include <complexes/NearestNeighbours.hh>

namespace aleph
{

namespace complexes
{

template <class Container> class FLANN : public NearestNeighbours<FLANN>
{
public:
  FLANN( const Container& container )
  {
  }

private:
  const Container& _container;
};

}

}

#endif
