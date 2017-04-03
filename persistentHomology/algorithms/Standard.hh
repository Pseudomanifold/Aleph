#ifndef ALEPH_PERSISTENT_HOMOLOGY_ALGORITHMS_STANDARD_HH__
#define ALEPH_PERSISTENT_HOMOLOGY_ALGORITHMS_STANDARD_HH__

#include "topology/BoundaryMatrix.hh"

#include <tuple>
#include <vector>

namespace aleph
{

namespace persistentHomology
{

namespace algorithms
{

class Standard
{
public:
  template <class Representation> void operator()( topology::BoundaryMatrix<Representation>& M )
  {
    using Index = typename Representation::Index;

    auto numColumns = M.getNumColumns();

    std::vector< std::pair<Index, bool> > lut( static_cast<std::size_t>( numColumns ),
                                               std::make_pair(0, false) );

    for( Index j = 0; j < numColumns; j++ )
    {
      Index i;
      bool valid = false;

      std::tie( i, valid ) = M.getMaximumIndex( j );
      while( valid && lut[ static_cast<std::size_t>(i) ].second )
      {
        M.addColumns( lut[ static_cast<std::size_t>(i) ].first, j );
        std::tie( i, valid ) = M.getMaximumIndex( j );
      }

      if( valid )
        lut[ static_cast<std::size_t>(i) ] = std::make_pair( j, true );
    }
  }
};

} // namespace algorithms

} // namespace persistentHomology

} // namespace aleph

#endif
