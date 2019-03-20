#ifndef ALEPH_TOPOLOGY_IO_ADJACENCY_MATRIX_HH__
#define ALEPH_TOPOLOGY_IO_ADJACENCY_MATRIX_HH__

#include <aleph/utilities/String.hh>

#include <aleph/topology/filtrations/Data.hh>

#include <algorithm>
#include <fstream>
#include <istream>
#include <string>
#include <unordered_map>
#include <vector>

#include <cmath>

namespace aleph
{

namespace topology
{

namespace io
{

/**
  @class AdjacencyMatrixReader
  @brief Reads square adjacency matrices in text format

  This reader class is meant to load square adjacency matrices (square)
  in text format. Entry (i,j) in the matrix contains the edge weight of
  the (unique) edge connecting nodes i and j.

  Depending on the configuration of the class, cells with a pre-defined
  weight (usually zero) are taken to indicate missing edges.

  The number of rows and columns must not vary over the file. An *empty*
  line is permitted, though. Likewise, lines starting with `#` will just
  be ignored. An example of a 3-by-3 matrix follows:

  \code
  0 1 2
  3 4 5
  2 1 7
  \endcode

  All simplicial complexes created by this class will be reported
  in filtration order, following the detected weights.
*/

class AdjacencyMatrixReader
{
public:

  /**
    Reads a simplicial complex from a file.

    @param filename Input filename
    @param  K       Simplicial complex
  */

  template <class SimplicialComplex> void operator()( const std::string& filename, SimplicialComplex& K )
  {
    std::ifstream in( filename );
    if( !in )
      throw std::runtime_error( "Unable to read input file" );

    this->operator()( in, K );
  }

  /** @overload operator()( const std::string&, SimplicialComplex& ) */
  template <class SimplicialComplex> void operator()( std::istream& in, SimplicialComplex& K )
  {
    using Simplex    = typename SimplicialComplex::ValueType;
    using DataType   = typename Simplex::DataType;
    using VertexType = typename Simplex::VertexType;

    auto position    = in.tellg();
    std::size_t n    = 0;

    // An 'unrolled' version of all edge weights that can be read from
    // the file. They are supposed to correspond to a matrix with some
    // number of columns and some number of rows.
    std::vector<DataType> values;

    using namespace aleph::utilities;

    {
      std::string line;
      while( std::getline( in, line ) )
      {
        // Skip empty lines and comments as promised
        if( line.empty() || line.front() == '#' )
          continue;

        ++n;
      }

      in.clear();
      in.seekg( position );
    }

    // FIXME: this will not work in case comments are part of the file.
    // Drat---should probably rewrite it.
    std::copy( std::istream_iterator<DataType>( in ), std::istream_iterator<DataType>(),
               std::back_inserter( values ) );

    // We cannot fill an empty simplicial complex. It might be useful to
    // throw an error here, though.
    if( values.empty() )
      return;

    _dimension = n;

    if( values.size() != _dimension * _dimension )
      throw std::runtime_error( "Format error: number of columns must not vary" );

    std::vector<Simplex> simplices;

    // First part of that equation reserves $n$ nodes, where $n$ is the
    // dimension of the matrix, followed by at most $n^2$ edges. Notice
    // that this assumes that *all* edges are defined.
    simplices.reserve( _dimension + ( _dimension * _dimension ) );

    // Edges -----------------------------------------------------------
    //
    // Create the edges first and update information about their weights
    // along with them.

    for( std::size_t y = 0; y < _dimension; y++ )
    {
      for( std::size_t x = 0; x < _dimension; x++ )
      {
        auto i = static_cast<VertexType>( _dimension * y + x   );
        auto w = values[i];

        // Map matrix indices to the corresponding vertex indices as
        // outlined above.
        auto u = VertexType(y);
        auto v = VertexType(x + _dimension);

        // TODO: skip edge if its weight is going to be ignored anyway
        // according to the client input.

        // We have no choice here but to store the corresponding simplex
        // with *exactly* the weight as it was specified in the file.
        simplices.push_back( Simplex( {u,v}, w ) );
      }
    }

    // Vertices --------------------------------------------------------
    //
    // Create a vertex for every node in the input data. This will use
    // the minimum weight detected in the file.
    auto minWeight = *std::min_element( values.begin(), values.end() );

    for( std::size_t i = 0; i < _dimension; i++ )
    {
      simplices.push_back(
        Simplex( VertexType( i ), minWeight )
      );
    }

    K = SimplicialComplex( simplices.begin(), simplices.end() );

    // Establish filtration order based on weights. There does not seem
    // to be much of a point to make this configurable; the edge weight
    // is a given property of the data.
    K.sort(
      filtrations::Data<Simplex>()
    );
  }

  /** @returns Dimension of matrix that was read last */
  std::size_t dimension() const noexcept { return _dimension; }

private:

  // Dimension of the matrix that was read last by this reader; this
  // will only be set if the matrix is actually square.
  std::size_t _dimension = 0;
};

} // namespace io

} // namespace topology

} // namespace aleph

#endif
