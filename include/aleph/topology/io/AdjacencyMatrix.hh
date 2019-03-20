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

    _height = height;
    _width  = values.size() / height;

    if( values.size() % height != 0 )
      throw std::runtime_error( "Format error: number of columns must not vary" );

    // This is required in order to assign the weight of nodes
    // correctly; we cannot trust the weights to be positive.
    auto minData = *std::min_element( values.begin(), values.end() );

    std::vector<Simplex> simplices;
    simplices.reserve( _width * _height + ( _width + _height ) );

    // Edges -----------------------------------------------------------
    //
    // Create the edges first and update information about their weights
    // along with them.

    // For determining the minimum weight, we first loop over all
    // possible edges, create a lookup table for the weights, and
    // finally create all the vertices using this lookup table. A
    // vertex will only information stored here if the right flag
    // has been set by the client.
    std::unordered_map<VertexType, DataType> minWeight;

    auto updateOrSetWeight
      = [&minWeight, this] ( const VertexType& v, const DataType& w )
        {
          if( minWeight.find( v ) == minWeight.end() )
            minWeight[v] = w;
          else
          {
            if( _assignMinimumAbsoluteVertexWeight )
            {
              if( std::abs( w ) < std::abs( minWeight[v] ) )
                minWeight[v] = w;
            }
            else
              minWeight[v] = std::min( minWeight[v], w );
          }
        };

    for( std::size_t y = 0; y < _height; y++ )
    {
      for( std::size_t x = 0; x < _width; x++ )
      {
        auto i = static_cast<VertexType>( _width * y + x   );
        auto w = values[i];

        // Map matrix indices to the corresponding vertex indices as
        // outline above.
        auto u = VertexType(y);
        auto v = VertexType(x + _height);

        updateOrSetWeight( u, w );
        updateOrSetWeight( v, w );

        simplices.push_back( Simplex( {u,v}, w ) );
      }
    }

    // Vertices --------------------------------------------------------
    //
    // Create a vertex for every node in the input data. An (n,m)-matrix
    // thus gives rise to n+m nodes.

    for( std::size_t i = 0; i < _height + _width; i++ )
    {
      // Notice that that `minWeight` map is guaranteed to contain the
      // weight (potentially signed) that corresponds to the vertex. A
      // different way of setting up the map depends on the flags that
      // are set by the client.
      simplices.push_back(
        Simplex( VertexType( i ),
          _assignMinimumVertexWeight || _assignMinimumAbsoluteVertexWeight
            ?  minWeight[ VertexType(i) ]
            :  minData )
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
