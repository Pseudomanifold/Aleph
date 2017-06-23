#ifndef ALEPH_TOPOLOGY_IO_EDGE_LISTS_HH__
#define ALEPH_TOPOLOGY_IO_EDGE_LISTS_HH__

#include <algorithm>
#include <fstream>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

#include <aleph/utilities/String.hh>

namespace aleph
{

namespace topology
{

namespace io
{

class EdgeListReader
{
public:

  bool readWeights() const noexcept { return _readWeights; }
  bool trimLines()   const noexcept { return _trimLines;   }

  void setReadWeights( bool value = true ) noexcept { _readWeights = value; }
  void setTrimLines( bool value = true )   noexcept { _trimLines = value; }

  template <class SimplicialComplex> void operator()( const std::string& filename, SimplicialComplex& K )
  {
    std::ifstream in( filename );
    if( !in )
      throw std::runtime_error( "Unable to read input file" );

    this->operator()( in, K );
  }

  template <class SimplicialComplex> void operator()( std::ifstream& in, SimplicialComplex& K )
  {
    using namespace utilities;

    using Simplex           = typename SimplicialComplex::ValueType;
    using DataType          = typename Simplex::DataType;
    using VertexType        = typename Simplex::VertexType;

    std::string line;

    std::set<Simplex> vertices;
    std::vector<Simplex> edges;

    while( in )
    {
      std::getline( in, line );

      if( _trimLines )
        line = trim( line );

      // TODO: Make this configurable and permit splitting by different
      // tokens such as commas
      auto tokens = split( line );

      // Skip empty lines and comments
      if( line.empty() || std::find( _commentTokens.begin(), _commentTokens.end(), line.front() ) != _commentTokens.end() )
        continue;

      if( tokens.size() >= 2 )
      {
        // TODO: Make order of vertices & weights configurable?
        VertexType u = convert<VertexType>( tokens[0] );
        VertexType v = convert<VertexType>( tokens[1] );
        DataType   w = DataType();

        if( tokens.size() >= 3 && _readWeights )
          w = convert<DataType>( tokens[2] );

        edges.push_back( Simplex( { u, v }, w ) );

        vertices.insert( Simplex( u ) );
        vertices.insert( Simplex( v ) );
      }
      else
      {
        // TODO: Throw error?
      }
    }

    // Using a set has the advantage of ensuring that duplicate
    // simplices are deleted automatically. A duplicate simplex
    // is usually created by the input data set itself. It must
    // not be considered for any subsequent analysis.
    std::set<Simplex> simplices;

    simplices.insert( vertices.begin(), vertices.end() );
    simplices.insert( edges.begin(), edges.end() );

    K = SimplicialComplex( simplices.begin(), simplices.end() );
  }

private:
  std::vector<char> _commentTokens = { '#', '%', '\"', '*' };

  bool _readWeights              = true;
  bool _trimLines                = true;
};

} // namespace io

} // namespace topology

} // namespace aleph

#endif
