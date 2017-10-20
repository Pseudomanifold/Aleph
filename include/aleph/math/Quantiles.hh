#ifndef ALEPH_MATH_QUANTILES_HH__
#define ALEPH_MATH_QUANTILES_HH__

#include <algorithm>
#include <iterator>

namespace aleph
{

namespace math
{

/** Calculates the median value of an input iterator range. */
template <class InputIterator> typename std::iterator_traits<InputIterator>::value_type median( InputIterator begin, InputIterator end )
{
  auto n = std::distance( begin, end );
  auto m = n / 2;

  // In order to avoid changing the original values, let's make a copy
  // of the range.
  using T = typename std::iterator_traits<InputIterator>::value_Type;
  std::vector<T> values( begin, end );

  std::nth_element( values.begin(), values.begin() + m, values.end() );

  auto middle = values[m];

  if( n % 2 )
    return middle;
  else
  {
    std::nth_element( values.begin(), values.begin() + m - 1, values.end() );
    return ( middle + values[m-1] ) / T(2);
  }
}

} // namespace math

} // namespace aleph

#endif
