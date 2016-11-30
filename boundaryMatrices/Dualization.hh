#ifndef ALEPH_DUALIZATION_HH__
#define ALEPH_DUALIZATION_HH__

#include "BoundaryMatrix.hh"

#include <algorithm>

namespace aleph
{

template <class Representation> BoundaryMatrix<Representation> dualize( const BoundaryMatrix<Representation>& M )
{
  using Index = typename BoundaryMatrix<Representation>::Index;

  auto&& numColumns = M.getNumColumns();

  std::vector< std::vector<Index> > dualMatrix( numColumns );
  std::vector<Index> dualDimensions( numColumns );
  std::vector<std::size_t> dualColumnSizes( numColumns );

  // Determine the size of every column in the dualized matrix. This
  // keeps memory re-allocation at a minimum.

  for( Index j = 0; j < numColumns; j++ )
  {
    auto&& column = M.getColumn(j);

    for( auto&& i : column )
      ++dualColumnSizes[ numColumns - 1 - i ];
  }

  for( Index j = 0; j < numColumns; j++ )
    dualMatrix[j].reserve( dualColumnSizes[j] );

  // Calculate the actual anti-transpose of the matrix. Since the
  // vectors have been properly resized, this operation should be
  // relatively well-behaved.

  for( Index j = 0; j < numColumns; j++ )
  {
    auto&& column = M.getColumn( j );

    for( auto&& i : column )
      dualMatrix[ numColumns - 1 - i ].push_back( numColumns - 1 - j );
  }

  auto&& d = M.getDimension();

  // FIXME: Do I need this?
  for( Index j = 0; j < numColumns; j++ )
    dualDimensions.at( numColumns - 1 - j ) = d - M.getDimension( j ); // FIXME: Change operator later after debugging

  BoundaryMatrix<Representation> N;
  N.setNumColumns( static_cast<Index>( dualMatrix.size() ) );

  for( Index j = 0; j < N.getNumColumns(); j++ )
  {
    // Do not assume that the column is properly sorted. A normal
    // std::reverse should be sufficient in most cases here but I
    // do not want to take any chances.
    std::sort( dualMatrix.at(j).begin(), dualMatrix.at(j).end() );

    N.setColumn( j,
                 dualMatrix.at(j).begin(), dualMatrix.at(j).end() );
  }

  N.setDualized( !M.isDualized() );
  return N;
}

}

#endif
