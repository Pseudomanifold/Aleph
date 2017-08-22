#ifndef ALEPH_TOPOLOGY_SKELETON_HH__
#define ALEPH_TOPOLOGY_SKELETON_HH__

namespace aleph
{

namespace topology
{

/**
  Stateless functor for extracting the \f$k\f$-skeleton of a simplicial
  complex. The functor is supposed to be used like a function:

  \code{.cpp}
  aleph::topology::Skeleton skeleton;

  auto K0 = skeleton(0, K); // 0-skeleton
  auto K1 = skeleton(1, K); // 1-skeleton
  \endcode

  The order of the extraction respects the filtration order of the
  given simplicial complex.
*/

class Skeleton
{
public:

  /**
    Extracts the \f$k\f$-skeleton of a given simplicial complex. While
    the function maintains the initial order of the input complex, the
    returned skeleton may be empty.
  */

  template <class SimplicialComplex> SimplicialComplex operator()( std::size_t k, const SimplicialComplex& K ) const noexcept
  {
    SimplicialComplex L;

    for( auto&& simplex : K )
    {
      if( simplex.dimension() <= k )
        L.push_back( simplex );
    }

    return L;
  }
};

} // namespace topology

} // namespace aleph

#endif
