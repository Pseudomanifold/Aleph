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

  std::vector<Index> indices( 2*functionValues.size() - 1 );
  std::iota( indices.begin(), indices.end(), Index(0) );

  std::stable_sort( indices.begin(), indices.end(),
             [&functionValues] ( Index i, Index j )
             {
               auto weight = [&] ( Index i )
               {
                 if( i < functionValues.size() )
                   return functionValues.at(i);
                 else
                 {
                   auto l = functionValues.at( i - functionValues.size()     );
                   auto r = functionValues.at( i - functionValues.size() + 1 );

                   return std::max( l, r );
                 }
               };

               return weight(i) < weight(j);
             } );

  for( Index j = 0; j < static_cast<Index>( indices.size() ); j++ )
  {
    auto&& index = indices.at(j);

    if( index < functionValues.size() )
      boundaryMatrix.clearColumn( j );
    else
    {
      std::vector<Index> vertexIndices = { index, index + 1 };

      boundaryMatrix.setColumn(j,
                               vertexIndices.begin(), vertexIndices.end() );
    }
  }
}

}

}

#endif
