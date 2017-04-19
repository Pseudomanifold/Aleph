#ifndef ALEPH_MATH_KAHAN_SUMMATION_HH__
#define ALEPH_MATH_KAHAN_SUMMATION_HH__

#include <algorithm>
#include <vector>

namespace aleph
{

namespace math
{

template <class T> class KahanSummation
{
public:
  KahanSummation( T initial = T() )
    : _sum( initial )
    , _c( T() )
  {
  }

  KahanSummation& operator+= ( T v )
  {
    T y  = v - _c;
    T t  = _sum + y;
    _c   = ( t - _sum ) - y;
    _sum = t;

    return *this;
  }

  KahanSummation& operator-= ( T v )
  {
    return this->operator+=( -v );
  }

  KahanSummation& operator*= ( T v )
  {
    _sum *= v;
    return *this;
  }

  KahanSummation& operator/= ( T v )
  {
    _sum /= v;
    return *this;
  }

  operator const T () const
  {
    return _sum;
  }

private:
  T _sum;
  T _c;
};

// ---------------------------------------------------------------------

/**
  Accumulation function modelled after \c std::accumulate. Instead of summing
  up the values without regards to floating point cancellation, this function
  uses the Kahan summation algorithm on a sorted range of numbers. This gives
  better stability guarantees for the sum.
*/

template <class InputIterator, class T> T accumulate_kahan( InputIterator first, InputIterator last, T init )
{
  KahanSummation<T> sum = init;

  std::vector<T> values( first, last );
  std::sort( values.begin(), values.end() );

  for( auto&& value : values )
    sum += value;

  return sum;
}

} // namespace math

} // namespace aleph

#endif
