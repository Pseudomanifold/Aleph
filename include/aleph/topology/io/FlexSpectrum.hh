#ifndef ALEPH_TOPOLOGY_IO_FLEX_SPECTRUM_HH__
#define ALEPH_TOPOLOGY_IO_FLEX_SPECTRUM_HH__

#include <aleph/topology/filtrations/Data.hh>

#include <algorithm>
#include <fstream>
#include <map>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

#include <cctype>

namespace aleph
{

namespace topology
{

namespace io
{

class FlexSpectrumReader
{
public:
  template <class SimplicialComplex> void operator()( const std::string& filename, SimplicialComplex& K )
  {
    std::ifstream in( filename );
    if( !in )
      throw std::runtime_error( "Unable to read input file" );

    using Simplex    = typename SimplicialComplex::ValueType;
    using DataType   = typename Simplex::DataType;
    using VertexType = typename Simplex::VertexType;

    std::string line;
    std::size_t index = 0; // line index, but conveniently, this will
                           // also be the vertex index

    std::vector<DataType> intensities;
    std::vector<Simplex> simplices;

    // Read lines & create vertices ------------------------------------
    //
    // This only fills up the auxiliary data structure of the spectrum
    // and creates a set of vertices for the simplicial complex.

    while( std::getline( in, line ) )
    {
      // Skip all lines that we cannot read. The data format is assumed
      // to contain a digit at the beginning of the line. Otherwise, we
      // do not consider this to be a valid line that contains data.
      if( line.empty() || !std::isdigit( line.front() ) )
        continue;

      std::stringstream converter( line );

      double x   = 0.0;
      DataType y = DataType();

      converter >> x
                >> y;

      _index_to_value[index] = x;

      simplices.emplace_back( Simplex( VertexType( index ), y ) );
      intensities.emplace_back( y );

      ++index;
    }

    // Abandon further calculations; the simplicial complex is empty or
    // almost empty and no edges can be added.
    if( index <= 1)
    {
      K = SimplicialComplex( simplices.begin(), simplices.end() );
      return;
    }

    // Create edges ----------------------------------------------------
    //
    // Creates a superlevel set filtration for the spectrum by assigning
    // each edge the *minimum* of the two function values.

    for( std::size_t i = 0; i < index - 1; i++ )
    {
      VertexType u = VertexType(i  );
      VertexType v = VertexType(i+1);
      DataType   w = std::min( simplices.at(u).data(), simplices.at(v).data() );

      simplices.emplace_back( Simplex( {u,v}, w ) );
    }

    if( _normalize )
    {
      auto totalIntensity
        = std::accumulate( intensities.begin(), intensities.end(), DataType() );

      for( auto&& s : simplices )
        s.setData( s.data() / totalIntensity );
    }

    K = SimplicialComplex( simplices.begin(), simplices.end() );
    K.sort(
      aleph::topology::filtrations::Data<Simplex,
                                         std::greater<DataType> >()
    );
  }

  bool normalize() const noexcept
  {
    return _normalize;
  }

  void normalize( bool value = true )
  {
    _normalize = value;
  }

private:

  /**
    If set, normalizes each spectrum according to all intensities that
    have been observed. This ensures that all masses sum up to one.
  */

  bool _normalize = false;

  /**
    Contains the raw \f$x\f-values that have been read from the file
    during the last parsing process. This permits the reconstruction
    of the original data.
  */

  std::map<std::size_t, double> _index_to_value;
};

} // namespace io

} // namespace topology

} // namespace aleph

#endif
