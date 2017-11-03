#ifndef ALEPH_TOPOLOGY_FILTER_HH__
#define ALEPH_TOPOLOGY_FILTER_HH__

namespace aleph
{

namespace topology
{

/**
  @class Filter
  @brief Stateless filter class functor for simplicial complexes

  This class permits *filtering* a given simplicial complex using an
  arbitrary functor that acts as a predicate. A good example of such
  a filter is the *removal* of simplices according to some condition
  such as their dimension, weight, and so on.

  The functor has to satisfy the following interface:

  \code{.cpp}
  template <class Simplex> bool operator()( const Simplex& s )
  {
    // true signifies that the simplex is retained in the resulting
    // simplicial complex
    return true;
  }
  \endcode
*/

class Filter
{
public:
  template <class SimplicialComplex, class Functor> SimplicialComplex operator()( const SimplicialComplex& K, Functor f ) const noexcept
  {
    SimplicialComplex L;

    for( auto&& simplex : K )
    {
      if( f( simplex ) )
        L.push_back( simplex );
    }

    return L;
  }
};

} // namespace topology

} // namespace aleph

#endif
