#ifndef ALEPH_UTILITIES_EMPTY_FUNCTOR_HH__
#define ALEPH_UTILITIES_EMPTY_FUNCTOR_HH__

namespace aleph
{

namespace utilities
{

/**
  @class EmptyFunctor
  @brief Basic drop-in empty functor class

  This class is a generic drop-in replacement for functors taking an arbitrary
  number of arguments. This is useful when algorithms expect a functor but the
  client does not require providing one.
*/

class EmptyFunctor
{
public:
  template <class... T> bool operator()( T... /* arguments */ )
  {
    return true;
  }
};

} // namespace utilities

} // namespace aleph

#endif
