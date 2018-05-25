#ifndef ALEPH_UTILITIES_VALUES_HH__
#define ALEPH_UTILITIES_VALUES_HH__

#include <stdexcept>

namespace aleph
{

namespace utilities
{

/**
  Checks whether a value is within a given range and throws an exception
  if this is not the case.

  @param x     Value
  @param lower Lower bound that the value must satisfy
  @param upper Upper bound that the value must satisfy
*/

template <class T> void ensureRange( T x, T lower, T upper )
{
  bool valid = x >= lower && x <= upper;
  if( !valid )
    throw std::out_of_range( "Value is out of bounds" );
}

/**
  Checks whether a value is larger than a given bound and throws an
  exception if this is not the case.

  @param x     Value
  @param lower Lower bound that the value must satisfy
*/

template <class T> void ensureLarger( T x, T lower )
{
  bool valid = x > lower;
  if( !valid )
    throw std::out_of_range( "Value is out of bounds" );
}

} // namespace utilities

} // namespace aleph

#endif
