#ifndef ALEPH_UTILITIES_UNORDERED_SET_OPERATIONS_HH__
#define ALEPH_UTILITIES_UNORDERED_SET_OPERATIONS_HH__

#include <unordered_set>

namespace aleph
{

namespace utilities
{

template <class T> void set_intersection( const std::unordered_set<T>& A,
                                          const std::unordered_set<T>& B,
                                          std::unordered_set<T>& C )
{
  // Let's ensure that we always use the smaller set in order to iterate over
  // the elements in the set.
  if( B.size() < A.size() )
  {
    set_intersection( B, A, C );
    return;
  }

  for( auto&& a : A )
  {
    if( B.find(a) != B.end() )
      C.insert(a);
  }
}

} // namespace utilities

} // namespace aleph

#endif
