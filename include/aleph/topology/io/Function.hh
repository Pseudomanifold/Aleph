#ifndef ALEPH_TOPOLOGY_IO_FUNCTION_HH__
#define ALEPH_TOPOLOGY_IO_FUNCTION_HH__

#include <aleph/topology/BoundaryMatrix.hh>

#include <algorithm>
#include <fstream>
#include <iterator>
#include <numeric>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include <iostream>

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

/**
  Converts function values to a simplicial complex. The client
  can define the corresponding vertex and data types. Functors
  can be used to determine how weights are being calculated.
*/

template <
  class SimplicialComplex,
  class Functor,
  class InputIterator
> SimplicialComplex loadFunction( InputIterator begin, InputIterator end,
                                  Functor f )
{
  using Simplex    = typename SimplicialComplex::ValueType;
  using DataType   = typename Simplex::DataType;
  using VertexType = typename Simplex::VertexType;

  static_assert( std::is_same<DataType, decltype( f(DataType(), DataType() ) )>::value,
                 "Functor return type must be compatible with simplex data type" );

  SimplicialComplex K;

  VertexType vertex = VertexType();
  for( auto it = begin; it != end; ++it )
    K.push_back( Simplex( vertex++, *it ) );

  vertex = VertexType();
  for( auto it = begin; it != end; ++it )
  {
    auto next = std::next( it );
    if( next == end )
      break;

    auto&& a = *it;
    auto&& b = *next;
    auto   w = f(a,b);

    // This is an edge that connects two adjacent simplices; the
    // weight is set according to their maximum by default.
    K.push_back( Simplex( {vertex, VertexType( vertex+1 ) }, w ) );

    ++vertex;
  }

  return K;
}

/**
  Loads a set of 1D functions from a file and converts them to
  simplicial complexes. The file format is simple and consists
  of a number of function values per line. A line break starts
  a new function. For example:

  \verbatim
  0 1 2 3
  3 1 6 4
  \endverbatim

  The preceding block describes two functions with 4 vertices.
*/

template <class SimplicialComplex, class Functor> std::vector<SimplicialComplex> loadFunctions( const std::string& filename,
                                                                                                Functor f )
{
  std::ifstream in( filename );
  if( !in )
    throw std::runtime_error( "Unable to read input file" );

  using Simplex    = typename SimplicialComplex::ValueType;
  using DataType   = typename Simplex::DataType;
  using VertexType = typename Simplex::VertexType;

  static_assert( std::is_same<DataType, decltype( f(DataType(), DataType() ) )>::value,
                 "Functor return type must be compatible with simplex data type" );

  std::vector<SimplicialComplex> complexes;

  std::string line;
  while( std::getline( in, line ) )
  {
    std::vector<DataType> functionValues;

    std::istringstream stream( line );

    std::copy( std::istream_iterator<DataType>( stream ),
               std::istream_iterator<DataType>(),
               std::back_inserter( functionValues ) );

    SimplicialComplex K;

    VertexType vertex = VertexType();

    for( auto&& value : functionValues )
      K.push_back( Simplex( vertex++, value ) );

    vertex = VertexType();
    for( auto it = functionValues.begin(); it != functionValues.end(); ++it )
    {
      auto next = std::next( it );
      if( next == functionValues.end() )
        break;

      auto&& a = *it;
      auto&& b = *next;
      auto   w = f(a,b);

      // This is an edge that connects two adjacent simplices; the
      // weight is set according to their maximum by default.
      K.push_back( Simplex( {vertex, VertexType( vertex+1 ) }, w ) );

      ++vertex;
    }

    complexes.push_back( K );
  }

  return complexes;
}

template <class SimplicialComplex> std::vector<SimplicialComplex> loadFunctions( const std::string& filename )
{
  using Simplex    = typename SimplicialComplex::ValueType;
  using DataType   = typename Simplex::DataType;

  return loadFunctions<SimplicialComplex>( filename, [] ( DataType a, DataType b ) { return std::max(a,b); } );
}

} // namespace io

} // namespace topology

} // namespace aleph

#endif
