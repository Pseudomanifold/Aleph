#ifndef ALEPH_DISTANCES_TRAITS_HH__
#define ALEPH_DISTANCES_TRAITS_HH__

namespace aleph
{

namespace distances
{

/**
  A generic traits class for distance functors. The basic idea behind this class
  is that different distance functors, in particular those that are based on the
  $L_p$ distances, may internally use squared distances for easier calculations.

  When giving clients the option to use distance-based methods, though, they are
  expecting the distances to be given in an unmodified form. Hence, every traits
  class provides two methods for converting distances from and to the class.
*/
  
template <class T> struct Traits
{
  using ResultType  = typename T::ResultType;
  using ElementType = typename T::ElementType;

  ResultType from( ElementType x ) const noexcept
  {
    return ResultType( x );
  }

  ResultType to( ElementType x ) const noexcept
  {
    return ResultType( x );
  }
};

}

}

#endif
