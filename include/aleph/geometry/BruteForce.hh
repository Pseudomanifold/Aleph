#ifndef ALEPH_GEOMETRY_BRUTE_FORCE_HH__
#define ALEPH_GEOMETRY_BRUTE_FORCE_HH__

#include <aleph/geometry/NearestNeighbours.hh>
#include <aleph/geometry/distances/Traits.hh>

#include <algorithm>
#include <vector>

namespace aleph
{

namespace geometry
{

/**
  @class BruteForce
  @brief Permits brute-force calculation of nearest neighbours

  This is a fallback solution for when there are no other options
  available for the calculation of nearest neighbours. This class
  enumerates all pairs of points in order to determine those that
  are within the specified radius of each other.
*/

template <class Container, class DistanceFunctor>
class BruteForce : public NearestNeighbours< BruteForce<Container, DistanceFunctor>, std::size_t, typename Container::ElementType >
{
public:
  using IndexType       = std::size_t;
  using ElementType     = typename Container::ElementType;
  using Traits          = aleph::distances::Traits<DistanceFunctor>;
  using Distance        = DistanceFunctor;

  explicit BruteForce( const Container& container )
    : _container( container )
  {
  }

  void radiusSearch( ElementType radius,
                     std::vector< std::vector<IndexType> >& indices,
                     std::vector< std::vector<ElementType> >& distances ) const
  {
    indices.clear();
    distances.clear();

    indices.resize( this->size() );
    distances.resize( this->size() );

    auto D               = _container.dimension();
    DistanceFunctor dist = DistanceFunctor();

    for( IndexType i = 0; i < this->size(); i++ )
    {
      // I am not making any assumptions about the distance functor
      // here. If it is not symmetric---and hence not a metric---we
      // really need to traverse all pairs.
      for( IndexType j = 0; j < this->size(); j++ )
      {
        auto d = dist( _container[i].begin(),
                       _container[j].begin(),
                       D );

        d = _traits.from( d );

        if( d < radius )
        {
          indices[i].push_back( j );
          distances[i].push_back( d );
        }
      }
    }
  }

  void neighbourSearch( unsigned k,
                        std::vector< std::vector<IndexType> >& indices,
                        std::vector< std::vector<ElementType> >& distances ) const
  {
    indices.clear();
    distances.clear();

    indices.resize( this->size() );
    distances.resize( this->size() );

    auto D               = _container.dimension();
    DistanceFunctor dist = DistanceFunctor();

    for( IndexType i = 0; i < this->size(); i++ )
    {
      // I am not making any assumptions about the distance functor
      // here. If it is not symmetric---and hence not a metric---we
      // really need to traverse all pairs.
      for( IndexType j = 0; j < this->size(); j++ )
      {
        auto d = dist( _container[i].begin(),
                       _container[j].begin(),
                       D );

        d = _traits.from( d );

        indices[i].push_back( j );
        distances[i].push_back( d );
      }

      std::sort( indices[i].begin(), indices[i].end(),
                 [&distances, &i] ( IndexType k, IndexType l )
                 {
                   return distances[i][k] < distances[i][l];
                 }
      );

      std::sort( distances[i].begin(), distances[i].end() );

      distances[i].erase( distances[i].begin() + k, distances[i].end() );
      indices[i].erase( indices[i].begin() + k, indices[i].end() );
    }
  }

  std::size_t size() const noexcept
  {
    return _container.size();
  }

private:

  /** Reference to the original container */
  const Container& _container;

  /** Required for optional distance functor conversions */
  Traits _traits;
};

} // namespace geometry

} // namespace aleph

#endif
