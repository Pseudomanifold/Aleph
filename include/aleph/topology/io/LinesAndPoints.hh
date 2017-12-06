#ifndef ALEPH_TOPOLOGY_IO_LINES_AND_POINTS_HH__
#define ALEPH_TOPOLOGY_IO_LINES_AND_POINTS_HH__

#include <fstream>
#include <ostream>
#include <stdexcept>
#include <string>

namespace aleph
{

namespace topology
{

namespace io
{

/**
  @class LinesAndPoints
  @brief Writes files as lines and points

  This writer class uses coordinate information of a simplicial complex
  to represent its 1-skeleton, i.e. its vertices and edges. All data is
  stored in a simple format, representing 0-simplex by its coordinates,
  and a 1-simplex by a sequence of coordinates.

  An edge \f$(u,v)\f$ will be stored as:

  \[
    x_u y_u
    x_v y_v\n
  \]

  Vertices and edges will be separated by two newlines in order to show
  that a new block was started. This makes it possible to use `gnuplot`
  for the visualization:

  \code
  plot "output.txt" index 0 with points pt 7, # to represent vertices
       ""           index 1 with lines        # to represent edges
  \endcode
*/

class LinesAndPoints
{
public:

  /** Writes the simplicial complex to a new file. */
  template <class SimplicialComplex, class Container> void operator()( const std::string& filename,
                                                                       const SimplicialComplex& K,
                                                                       const Container& container )
  {
    std::ofstream out( filename );
    if( !out )
      throw std::runtime_error( "Unable to open output file" );

    this->operator()( out, K, container );
  }

  /** Writes the simplicial complex to a new output stream. */
  template <class SimplicialComplex, class Container> void operator()( std::ostream& out,
                                                                       const SimplicialComplex& K,
                                                                       const Container& container )
  {
    for( auto&& s : K )
    {
      if( s.dimension() != 0 )
        continue;

      auto u = s[0];         // vertex
      auto p = container[u]; // coordinate

      if( p.size() >= 2 )
        write( out, p, _addVertexLabels ? std::to_string( u ) : std::string() );
      else
        throw std::runtime_error( "Insufficient number of dimensions for storing coordinates" );
    }

    out << "\n\n";

    for( auto&& s : K )
    {
      if( s.dimension() != 1 )
        continue;

      auto u = s[0];         // vertex
      auto v = s[1];         // vertex
      auto p = container[u]; // first coordinate
      auto q = container[v]; // second coordinate

      if( p.size() >= 2 && q.size() >= 2 )
      {
        write(out, p);
        write(out, q);

        out << "\n";
      }
      else
        throw std::runtime_error( "Insufficient number of dimensions for storing coordinates" );
    }
  }

  bool addVertexLabels() const noexcept
  {
    return _addVertexLabels;
  }

  void addVertexLabels( bool value ) noexcept
  {
    _addVertexLabels = value;
  }

private:
  template <class Container> void write( std::ostream& out, Container& p, const std::string& label = std::string() )
  {
    out << p[0] << "\t" << p[1];

    if( p.size() >= 3 )
      out << " " << p[2];

    if( not label.empty() )
      out << " " << label;

    out << "\n";
  }

  bool _addVertexLabels = false;
};

} // namespace io

} // namespace topology

} // namespace aleph


#endif
