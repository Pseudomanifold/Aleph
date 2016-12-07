#ifndef ALEPH_TOPOLOGY_IO_FUNCTION_HH__
#define ALEPH_TOPOLOGY_IO_FUNCTION_HH__

#include "topology/BoundaryMatrix.hh"

#include <algorithm>
#include <fstream>
#include <iterator>
#include <numeric>
#include <string>
#include <unordered_map>
#include <vector>

namespace aleph
{

namespace topology
{

namespace io
{

template <
  class DataType,
  class Representation
> void loadFunction( const std::string& filename,
                     BoundaryMatrix<Representation>& boundaryMatrix,
                     std::vector<DataType>& functionValues )
{
  using BoundaryMatrix = BoundaryMatrix<Representation>;
  using Index          = typename BoundaryMatrix::Index;

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

  boundaryMatrix.setNumColumns( static_cast<Index>( indices.size() ) );

  // Maps a vertex in the original function to its place in the current
  // filtration order. The map is filled while creating the matrix below,
  // which is possible because faces need to precede cofcaces.
  std::unordered_map<Index, Index> vertexIndexMap;

  for( Index j = 0; j < static_cast<Index>( indices.size() ); j++ )
  {
    auto&& index = indices.at(j);

    if( index < functionValues.size() )
    {
      boundaryMatrix.clearColumn( j );

      vertexIndexMap[index] = j;
    }
    else
    {
      Index k = static_cast<Index>( index - functionValues.size() );

      std::vector<Index> vertexIndices = { vertexIndexMap.at(k), vertexIndexMap.at(k+1) };

      boundaryMatrix.setColumn(j,
                               vertexIndices.begin(), vertexIndices.end() );
    }
  }

  // Extend function values with edge weights. Since the i-th edge will use
  // vertices i and i+1, this can be done in one sweep.
  {
    std::size_t n = functionValues.size();
    for( std::size_t i = 0; i < n - 1; i++ )
    {
      auto w1 = functionValues.at(i  );
      auto w2 = functionValues.at(i+1);
      functionValues.push_back( std::max( w1, w2 ) );
    }
  }

  // Sort the function values to reflect the order of vertex indices. Else, we
  // will be unable to add the proper weights to the corresponding pairing.

  {
    std::vector<DataType> newFunctionValues;
    newFunctionValues.reserve( functionValues.size() );

    for( auto&& index : indices )
      newFunctionValues.push_back( functionValues.at( index ) );

    functionValues.swap( newFunctionValues );
  }
}

} // namespace io

} // namespace topology

} // namespace aleph

#endif
