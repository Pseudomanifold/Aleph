#ifndef ALEPH_TOPOLOGY_IO_EDGE_LISTS_HH__
#define ALEPH_TOPOLOGY_IO_EDGE_LISTS_HH__

#include <algorithm>
#include <fstream>
#include <map>
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

/**
  @class EdgeListReader
  @brief Reader class for unstructured edge lists

  Graphs are often specified in this informal format in which every line
  consists of two vertex indices, separated by white-space, and followed
  by an optional weight. This reader permits loading these files in some
  variations. For example, weights may optionally not be read. Also, the
  *separator* that is used to separate vertex indices is configurable.

  Note that by default, this reader is capable of reading files in which
  different indices are separated by white-space. The indices may either
  be *numeric*, in which case they are automatically used, or *strings*,
  in which case they are converted to numerical IDs, following the order
  in which they occur in the file.

  The reader expects that vertex IDs are *followed* by weights per line,
  and considers `#`, `%`, `"`, `*` to be comment tokens while *skipping*
  empty lines.
*/

class EdgeListReader
{
public:

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

    std::size_t lastID = 0;

    while( in )
    {
      std::getline( in, line );

      if( _trimLines )
        line = trim( line );

      auto tokens = split( line, "[" + _separator + "]+" );

      // Skip empty lines and comments
      if( line.empty() || std::find( _commentTokens.begin(), _commentTokens.end(), line.front() ) != _commentTokens.end() )
        continue;

      if( tokens.size() >= 2 )
      {
        VertexType u = VertexType();
        VertexType v = VertexType();

        if( tokens[0].find_first_not_of( "0123456789" ) == std::string::npos )
        {
          u = convert<VertexType>( tokens[0] );
          v = convert<VertexType>( tokens[1] );
        }
        else
        {
          auto&& us = tokens[0];
          auto&& vs = tokens[1];

          if( _nodeLabels.find(us) == _nodeLabels.end() )
            _nodeLabels[us] = lastID++;

          if( _nodeLabels.find(vs) == _nodeLabels.end() )
            _nodeLabels[vs] = lastID++;

          u = static_cast<VertexType>( _nodeLabels.at(us) );
          v = static_cast<VertexType>( _nodeLabels.at(vs) );
        }

        DataType w = DataType();
        if( tokens.size() >= 3 && _readWeights )
          w = convert<DataType>( tokens[2] );

        edges.push_back( Simplex( { u, v }, w ) );

        vertices.insert( Simplex( u ) );
        vertices.insert( Simplex( v ) );
      }
      else
        throw std::runtime_error(" Format error: not enough tokens to continue parsing" );
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

  bool readWeights() const noexcept { return _readWeights; }
  bool trimLines()   const noexcept { return _trimLines;   }

  void setReadWeights( bool value = true ) noexcept { _readWeights = value; }
  void setTrimLines( bool value = true )   noexcept { _trimLines = value; }

  /**
    Sets the separator to use for splitting tokens on every line of an
    input file. Setting this to ":", for example, means that different
    vertex indices for an edge are separated by ":" instead of space.

    Set this to the special value `[:space:]` to permit splitting by
    any white-space character. This is the default behaviour.

    @param separator Separator to use
  */

  void setSeparator( const std::string& separator )
  {
    _separator = separator;
  }

private:
  std::vector<char> _commentTokens = { '#', '%', '\"', '*' };
  std::string _separator           = "[:space:]";

  /**
    Optional set of node labels. These are read automatically in case an
    input data set does not use numeric labels. The idea is to map an ID
    in the form of a string to an index in corresponding graph.
  */

  std::map<std::string, std::size_t> _nodeLabels;

  bool _readWeights              = true;
  bool _trimLines                = true;
};

} // namespace io

} // namespace topology

} // namespace aleph

#endif
