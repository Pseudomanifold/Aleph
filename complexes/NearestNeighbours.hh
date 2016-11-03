#ifndef ALEPH_COMPLEXES_NEAREST_NEIGHBOURS_HH__
#define ALEPH_COMPLEXES_NEAREST_NEIGHBOURS_HH__

namespace aleph
{

namespace complexes
{

template <class Wrapper, class Container> class NearestNeighbours
{
public:
  NearestNeighbours( const Container& container )
    : _wrapper( container )
  {
  }

private:
  Wrapper _wrapper;
};

}

}

#endif
