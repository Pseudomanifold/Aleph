#ifndef ALEPH_FUNCTION_HH__
#define ALEPH_FUNCTION_HH__

#include "BoundaryMatrix.hh"

#include <algorithm>
#include <fstream>
#include <iterator>
#include <numeric>
#include <string>
#include <vector>

namespace aleph
{

namespace io
{

template <
  class Index,
  class DataType
> void loadFunction( const std::string& filename,
                     BoundaryMatrix<Index>& boundaryMatrix,
                     std::vector<DataType>& functionValues )
{
  using BoundaryMatrix = BoundaryMatrix<Index>;

  std::ifstream in( filename );

  if( !in )
    throw std::runtime_error( "Unable to open input filename" );

  functionValues.clear();
  functionValues.shrink_to_fit();

  std::copy( std::istream_iterator<DataType>( in ),
             std::istream_iterator<DataType>(),
             std::back_inserter( functionValues ) );

  if( functionValues.empty() )
    throw std::runtime_error( "Unable to load any function values" );

  std::vector<Index> vertexIndices( functionValues.size() );

  std::iota( vertexIndices.begin(), vertexIndices.end(),
             Index(0) );

  std::sort( vertexIndices.begin(), vertexIndices.end(),
             [&functionValues] ( Index i, Index j )
             {
              return functionValues.at(i) < functionValues.at(j);
             } );

  std::vector<Index> edgeIndices( functionValues.size() - 1 );

  std::iota( edgeIndices.begin(),
             edgeIndices.end(),
             static_cast<Index>( functionValues.size() ) );

  std::sort( edgeIndices.begin(), edgeIndices.end(),
             [&functionValues] ( Index i, Index j )
             {
              auto k = i - functionValues.size();
              auto l = j - functionValues.size();

              return functionValues.at(k) < functionValues.at(l);
             } );

  boundaryMatrix.setNumColumns( 2 * functionValues.size() - 1 );

  for( Index j = 0; j < static_cast<Index>( functionValues.size() ); j++ )
    boundaryMatrix.clearColumn( j );
}

}

}

#endif
