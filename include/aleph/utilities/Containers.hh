#ifndef ALEPH_UTILITIES_CONTAINERS_HH__
#define ALEPH_UTILITIES_CONTAINERS_HH__

#include <vector>

#include <cassert>

namespace aleph
{

namespace utilities
{

namespace std
{

struct minus
{
  template <class T> std::vector<T> operator()( const std::vector<T>& lhs, const std::vector<T>& rhs ) const
  {
    assert( lhs.size() == rhs.size() );

    std::vector<T> result;
    result.reserve( lhs.size() );

    for( auto it1 = lhs.begin(), it2 = rhs.begin(); it1 != lhs.end() && it2 != rhs.end(); ++it1, ++it2 )
      result.push_back( *it1 - *it2 );

    return result;
  }
}

} // namespace std

} // namespace utilities

} // namespace aleph

#endif
