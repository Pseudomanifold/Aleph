#ifndef ALEPH_COMPLEXES_FLANN_HH__
#define ALEPH_COMPLEXES_FLANN_HH__

#include <complexes/NearestNeighbours.hh>

#include <flann/flann.hpp>

namespace aleph
{

namespace complexes
{

template <class Container> class FLANN : public NearestNeighbours< FLANN<Container> >
{
public:
  using ElementType     = typename Container::ElementType;
  using DistanceFunctor = flann::L2<ElementType>;

  FLANN( const Container& container )
    : _container( container )
  {
    _matrix
      = flann::Matrix<ElementType>( container.data(),
                                    container.size(), container.dimension() );

    flann::IndexParams indexParameters
      = flann::KDTreeSingleIndexParams();

    _index
      = new flann::Index<DistanceFunctor>( _matrix, indexParameters );
  }

  ~FLANN()
  {
    delete _index;
  }

private:
  const Container& _container;

  /**
    Copy of container data. This makes interfacing with FLANN easier, at
    the expense of having large storage costs.
  */

  flann::Matrix<ElementType> _matrix;

  /** Index structure for queries. TODO: Make configurable/generic. */
  flann::Index<DistanceFunctor>* _index = nullptr;
};

}

}

#endif
