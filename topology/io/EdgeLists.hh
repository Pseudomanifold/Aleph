#ifndef ALEPH_TOPOLOGY_IO_EDGE_LISTS_HH__
#define ALEPH_TOPOLOGY_IO_EDGE_LISTS_HH__

#include <algorithm>
#include <fstream>
#include <set>
#include <string>
#include <vector>

#include "utilities/String.hh"

namespace aleph
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

    std::vector<Simplex> simplices;
    simplices.reserve( vertices.size() + edges.size() );

    simplices.insert( simplices.end(), vertices.begin(), vertices.end() );
    simplices.insert( simplices.end(), edges.begin(), edges.end() );

    K = SimplicialComplex( simplices.begin(), simplices.end() );
  }

private:
  std::vector<char> _commentTokens = { '#', '%', '\"', '*' };

  bool _readWeights              = true;
  bool _trimLines                = true;
};

}

}

#endif
