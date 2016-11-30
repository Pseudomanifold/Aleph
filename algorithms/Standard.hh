#ifndef ALEPH_ALGORITHMS_STANDARD_HH__
#define ALEPH_ALGORITHMS_STANDARD_HH__

#include "boundaryMatrices/BoundaryMatrix.hh"

#include <tuple>
#include <vector>

namespace aleph
{

namespace algorithms
{

class Standard
{
public:
  template <class Representation> void operator()( BoundaryMatrix<Representation>& M )
  {
    using Index = typename Representation::Index;

    auto numColumns = M.getNumColumns();

    std::vector< std::pair<Index, bool> > lut( numColumns,
                                               std::make_pair(0, false) );

    for( Index j = 0; j < numColumns; j++ )
    {
      Index i;
      bool valid = false;

      std::tie( i, valid ) = M.getMaximumIndex( j );
      while( valid && lut[i].second )
      {
        M.addColumns( lut[i].first, j );
        std::tie( i, valid ) = M.getMaximumIndex( j );
      }

      if( valid )
        lut[i] = std::make_pair( j, true );
    }
  }
};

}

}

#endif
