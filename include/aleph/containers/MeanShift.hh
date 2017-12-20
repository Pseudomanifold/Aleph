#ifndef ALEPH_CONTAINERS_MEAN_SHIFT_HH__
#define ALEPH_CONTAINERS_MEAN_SHIFT_HH__

#include <algorithm>
#include <iterator>
#include <vector>

namespace aleph
{

namespace containers
{

template <
  class Wrapper      ,
  class Container    ,
  class InputIterator,
  class OutputIterator
> void meanShiftSmoothing(
  const Container& container,
  InputIterator begin, InputIterator end,
  OutputIterator result,
  unsigned k,
  unsigned n = 1
)
{
  using T = typename std::iterator_traits<InputIterator>::value_type;
  auto N  = container.size();

  // Makes it easier to permit random access to the container; I am
  // assuming that the indices correspond to each other.
  std::vector<T> data( begin, end );

  Wrapper nearestNeighbours( container );

  using IndexType   = typename Wrapper::IndexType;
  using ElementType = typename Wrapper::ElementType;

  std::vector< std::vector<IndexType> > indices;
  std::vector< std::vector<ElementType> > distances;

  nearestNeighbours.neighbourSearch( k+1, indices, distances );

  for( unsigned iteration = 0; iteration < n; iteration++ )
  {
    std::vector<T> data_( data.size() );

    for( std::size_t i = 0; i < N; i++ )
    {
      auto&& neighbours_  = indices[i];
      auto&& distances_   = distances[i];

      // TODO: make configurable; it might also make sense to permit the
      // usage of arbitrary functors here
      aleph::math::KahanSummation<T> value        = 0.0;
      aleph::math::KahanSummation<T> sumOfWeights = 0.0;

      for( std::size_t j = 0; j < neighbours_.size(); j++ )
      {
        auto index    = neighbours_[j];
        auto weight   = distances_[j] > 0 ? 1.0 / ( distances_[j] * distances_[j] ) : 1.0;
        value        += data[ index ] * weight; // use data values from *previous* step to
                                                // perform the smoothing!
        sumOfWeights += weight;
      }

      value    /= sumOfWeights;
      data_[i]  = value;
    }

    data.swap( data_ );
  }

  std::copy( data.begin(), data.end(), result );
}

} // namespace containers

} // namespace aleph

#endif
