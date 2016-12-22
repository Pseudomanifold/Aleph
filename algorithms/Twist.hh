#ifndef ALEPH_ALGORITHMS_TWIST_HH__
#define ALEPH_ALGORITHMS_TWIST_HH__

#include "topology/BoundaryMatrix.hh"

#include <tuple>
#include <vector>

namespace aleph
{

namespace algorithms
{

class Twist
{
public:
  template <class Representation> void operator()( topology::BoundaryMatrix<Representation>& M )
  {
    using Index = typename Representation::Index;

    auto dimension  = M.getDimension();
    auto numColumns = M.getNumColumns();

    std::vector< std::pair<Index, bool> > lut( numColumns,
                                               std::make_pair(0, false) );

    for( Index d = dimension; d >= 1; d-- )
    {
      for( Index j = 0; j < numColumns; j++ )
      {
        if( M.getDimension( j ) == d )
        {
          Index i;
          bool valid = false;

          std::tie( i, valid ) = M.getMaximumIndex( j );
          while( valid && lut[ std::size_t(i) ].second )
          {
            M.addColumns( lut[ std::size_t(i) ].first, j );
            std::tie( i, valid ) = M.getMaximumIndex( j );
          }

          if( valid )
          {
            lut[ std::size_t(i) ] = std::make_pair( j, true );
            M.clearColumn( i );
          }
        }
      }
    }
  }
};

}

}

#endif
