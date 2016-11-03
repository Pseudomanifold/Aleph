#ifndef ALEPH_COMPLEXES_FLANN_HH__
#define ALEPH_COMPLEXES_FLANN_HH__

#include <complexes/NearestNeighbours.hh>

#include <flann/flann.hpp>

#include <algorithm>
#include <vector>

namespace aleph
{

namespace complexes
{

template <class Container> class FLANN : public NearestNeighbours< FLANN<Container>, std::size_t, typename Container::ElementType >
{
public:
  using IndexType       = std::size_t;
  using ElementType     = typename Container::ElementType;

  // TODO: Make configurable...
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

  void radiusSearch( ElementType radius,
                     std::vector< std::vector<IndexType> >& indices,
                     std::vector< std::vector<ElementType> >& distances )
  {

    flann::SearchParams searchParameters = flann::SearchParams();
    searchParameters.checks = flann::FLANN_CHECKS_UNLIMITED;

    std::vector< std::vector<int> > internalIndices;

    _index->radiusSearch( _matrix,
                          internalIndices,
                          distances,
                          radius,
                          searchParameters );

    // Perform transformation of indices -------------------------------

    indices.clear();
    indices.resize( _matrix.rows );

    for( std::size_t i = 0; i < internalIndices.size(); i++ )
    {
      indices[i] = std::vector<IndexType>( internalIndices[i].size() );

      std::transform( internalIndices[i].begin(), internalIndices[i].end(),
                      indices[i].begin(),
                      [] ( int j )
                      {
                        return static_cast<IndexType>( j );
                      } );
    }
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
