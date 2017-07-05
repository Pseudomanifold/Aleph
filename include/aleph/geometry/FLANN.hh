#ifndef ALEPH_GEOMETRY_FLANN_HH__
#define ALEPH_GEOMETRY_FLANN_HH__

#include <aleph/config/FLANN.hh>

#include <aleph/geometry/NearestNeighbours.hh>
#include <aleph/geometry/distances/Traits.hh>

#ifdef ALEPH_WITH_FLANN
  #include <flann/flann.hpp>
#endif

#include <algorithm>
#include <vector>

namespace aleph
{

namespace geometry
{

template <class Container, class DistanceFunctor>
class FLANN : public NearestNeighbours< FLANN<Container, DistanceFunctor>, std::size_t, typename Container::ElementType >
{
public:
  using IndexType       = std::size_t;
  using ElementType     = typename Container::ElementType;
  using Traits          = aleph::distances::Traits<DistanceFunctor>;

  explicit FLANN( const Container& container )
    : _container( container )
  {
#ifdef ALEPH_WITH_FLANN
    _matrix
      = flann::Matrix<ElementType>( container.data(),
                                    container.size(), container.dimension() );

    flann::IndexParams indexParameters
      = flann::KDTreeSingleIndexParams();

    _index
      = new flann::Index<DistanceFunctor>( _matrix, indexParameters );

    _index->buildIndex();
#endif
  }

  ~FLANN()
  {
#ifdef ALEPH_WITH_FLANN
    delete _index;
#endif
  }

#ifdef ALEPH_WITH_FLANN
  void radiusSearch( ElementType radius,
                     std::vector< std::vector<IndexType> >& indices,
                     std::vector< std::vector<ElementType> >& distances ) const
  {
    flann::SearchParams searchParameters = flann::SearchParams();
    searchParameters.checks = flann::FLANN_CHECKS_UNLIMITED;

    using ResultType = typename DistanceFunctor::ResultType;

    std::vector< std::vector<int> > internalIndices;
    std::vector< std::vector<ResultType> > internalDistances;

    _index->radiusSearch( _matrix,
                          internalIndices,
                          internalDistances,
                          static_cast<float>( _traits.to( radius ) ),
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

    // Perform transformation of distances -----------------------------

    if( internalDistances.empty() )
      return;

    distances.clear();
    distances.resize( internalDistances.size(), std::vector<ElementType>( internalDistances.front().size() ) );

    {
      using size_type = typename std::vector< std::vector<ElementType> >::size_type;

      for( size_type row = 0; row < distances.size(); row++ )
        for( size_type col = 0; col < distances[row].size(); col++ )
          distances[row][col] = static_cast<ElementType>( internalDistances[row][col] );
    }

    for( auto&& D : distances )
    {
      std::transform( D.begin(), D.end(), D.begin(),
                      [this] ( ElementType x )
                      {
                        return _traits.from( x );
                      } );
    }
  }

#else
  void radiusSearch( ElementType /* radius */,
                     std::vector< std::vector<IndexType> >& /* indices */,
                     std::vector< std::vector<ElementType> >& /* distances */ ) const
  {
  }
#endif

#ifdef ALEPH_WITH_FLANN
  void neighbourSearch( unsigned k,
                        std::vector< std::vector<IndexType> >& indices,
                        std::vector< std::vector<ElementType> >& distances ) const
  {
    // FLANN does *not* like being queries for no neighbours at all, so
    // let's play nice.
    if( k == 0 )
    {
      indices.clear();
      distances.clear();

      indices.resize(   _matrix.rows );
      distances.resize( _matrix.rows );

      return;
    }

    flann::SearchParams searchParameters = flann::SearchParams();
    searchParameters.checks = flann::FLANN_CHECKS_UNLIMITED;

    using ResultType = typename DistanceFunctor::ResultType;

    std::vector< std::vector<int> > internalIndices;
    std::vector< std::vector<ResultType> > internalDistances;

    _index->knnSearch( _matrix,
                       internalIndices,
                       internalDistances,
                       k,
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

    // Perform transformation of distances -----------------------------

    if( internalDistances.empty() )
      return;

    distances.clear();
    distances.resize( internalDistances.size(), std::vector<ElementType>( internalDistances.front().size() ) );

    {
      using size_type = typename std::vector< std::vector<ElementType> >::size_type;

      for( size_type row = 0; row < distances.size(); row++ )
        for( size_type col = 0; col < distances[row].size(); col++ )
          distances[row][col] = static_cast<ElementType>( internalDistances[row][col] );
    }

    for( auto&& D : distances )
    {
      std::transform( D.begin(), D.end(), D.begin(),
                      [this] ( ElementType x )
                      {
                        return _traits.from( x );
                      } );
    }
  }

#else
  void neighbourSearch( unsigned /* k */,
                        std::vector< std::vector<IndexType> >& /* indices */,
                        std::vector< std::vector<ElementType> >& /* distances */ ) const
  {
  }
#endif

  std::size_t size() const noexcept
  {
    return _container.size();
  }

  // The wrapper must not be copied. Else, clients  will run afoul of memory
  // management issues.
  FLANN( const FLANN& other )            = delete;
  FLANN& operator=( const FLANN& other ) = delete;

private:
  const Container& _container;

#ifdef ALEPH_WITH_FLANN

  /**
    Copy of container data. This makes interfacing with FLANN easier, at
    the expense of having large storage costs.
  */

  flann::Matrix<ElementType> _matrix;

  /** Index structure for queries. */
  flann::Index<DistanceFunctor>* _index = nullptr;

#endif

  /** Required for optional distance functor conversions */
  Traits _traits;
};

} // namespace geometry

} // namespace aleph

#endif
