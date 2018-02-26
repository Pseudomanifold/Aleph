#ifndef ALEPH_TOPOLOGY_IO_MATRIX_HH__
#define ALEPH_TOPOLOGY_IO_MATRIX_HH__

#include <algorithm>
#include <fstream>
#include <istream>
#include <string>
#include <unordered_map>
#include <vector>

#include <aleph/utilities/String.hh>

#include <aleph/topology/filtrations/Data.hh>

namespace aleph
{

namespace topology
{

namespace io
{

/**
  @class BipartiteAdjacencyMatrixReader
  @brief Reads bipartite adjacency matrices in text format

  This reader class is meant to load bipartite adjacency matrices in
  text format. Every row of the matrix represents edges that connect
  nodes from the first class with nodes of the second class. Weights
  that are non-zero are used to indicate the presence of an edge.

  The number of rows and columns must not vary over the file. An *empty*
  line is permitted, though. Likewise, lines starting with `#` will just
  be ignored. An example of a 2-by-3 matrix follows:

  \code
  0 1 2
  3 4 5
  \endcode

  All simplicial complexes created by this class will be reported
  in filtration order, following the detected weights.
*/

class BipartiteAdjacencyMatrixReader
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

    auto position = in.tellg();

    std::size_t height = 0;
    std::size_t width  = 0;

    using namespace aleph::utilities;

    {
      std::string line;
      while( std::getline( in, line ) )
      {
        line        = trim( line );
        auto tokens = split( line );

        if( width == 0 )
          width = tokens.size();
        else if( width != tokens.size() )
          throw std::runtime_error( "Format error: number of columns must not vary" );

        ++height;
      }

      in.clear();
      in.seekg( position );
    }

    _height = height;
    _width  = width;

    using Simplex    = typename SimplicialComplex::ValueType;
    using DataType   = typename Simplex::DataType;
    using VertexType = typename Simplex::VertexType;

    std::vector<DataType> values;
    values.reserve( _height * _width );

    std::copy( std::istream_iterator<DataType>( in ), std::istream_iterator<DataType>(),
               std::back_inserter( values ) );

    // We cannot fill an empt simplicial complex. It might be useful to
    // throw an error here, though.
    if( values.empty() )
      return;

    // This is required in order to assign the weight of nodes
    // correctly; we cannot trust the weights to be positive.
    auto minData = *std::min_element( values.begin(), values.end() );

    std::vector<Simplex> simplices;

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
      = [&minWeight] ( const VertexType& v, const DataType& w )
        {
          if( minWeight.find( v ) == minWeight.end() )
            minWeight[v] = w;
          else
            minWeight[v] = std::min( minWeight[v], w );
        };

    for( std::size_t y = 0; y < _height; y++ )
    {
      for( std::size_t x = 0; x < _width; x++ )
      {
        auto i = static_cast<VertexType>( width * y + x   );
        auto w = values.at( i );

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
      simplices.push_back(
        Simplex( VertexType( i ),
          _assignMinimumVertexWeight
            ?  minWeight.at( VertexType(i) )
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

  /** @returns Height of matrix that was read last */
  std::size_t height() const noexcept { return _height; }

  /** @returns Width of matrix that was read last */
  std::size_t width()  const noexcept { return _width;  }

  /** Permits changing the behaviour of vertex weight assignment */
  void setAssignMinimumVertexWeight( bool value = true )
  {
    _assignMinimumVertexWeight = value;
  }

  /** @returns Flag whether vertex weights are assigned to be minimal or not */
  bool assignMinimumVertexWeight() const noexcept
  {
    return _assignMinimumVertexWeight;
  }

private:
  std::size_t _height = 0;
  std::size_t _width  = 0;

  /**
    If set, assigns the minimum vertex weight according to the minimum
    edge weight that is connected to the given vertex. Else, *all* the
    vertices will get a weight of zero.
  */

  bool _assignMinimumVertexWeight = false;
};

} // namespace io

} // namespace topology

} // namespace aleph

#endif
