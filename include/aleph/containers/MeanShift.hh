#ifndef ALEPH_CONTAINERS_MEAN_SHIFT_HH__
#define ALEPH_CONTAINERS_MEAN_SHIFT_HH__

#include <iterator>
#include <vector>

namespace aleph
{

namespace containers
{

template <
  class Container    ,
  class Wrapper      ,
  class InputIterator,
> void meanShiftSmoothing(
  const Container& container,
  unsigned k,
  InputIterator begin, InputIterator end )
{
  using T = typename std::iterator_traits<InputIterator>::value_type;
  auto n  = container.size();

  // Makes it easier to permit random access to the container.
  std::vector<T> data( begin, end );
}

} // namespace containers

} // namespace aleph

#endif
