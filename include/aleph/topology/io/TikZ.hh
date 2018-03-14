#ifndef ALEPH_TOPOLOGY_IO_TIKZ_HH__
#define ALEPH_TOPOLOGY_IO_TIKZ_HH__

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
  @class TikZ
  @brief Writes files as a TikZ picture

  This writer class uses coordinate information of a simplicial complex
  to represent its 1-skeleton, i.e. its vertices and edges. All data is
  stored as a picture for the TikZ LaTeX toolkit, resulting in a clean,
  modern output for publications.

  The output is configurable and permits the following adjustments:

    - point size
    - line width
*/

class TikZ
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
    out << "\\begin{tikzpicture}\n";

    // Render vertices as points ---------------------------------------

    out << "% 0-simplices\n";

    for( auto&& s : K )
    {
      if( s.dimension() != 0 )
        continue;

      auto u = s[0];         // vertex
      auto p = container[u]; // coordinate

      if( p.size() >= 2 )
        writePoint( out, u, p );
      else
        throw std::runtime_error( "Insufficient number of dimensions for storing coordinates" );
    }

    // Render edges as lines -------------------------------------------

    out << "% 1-simplices\n";

    for( auto&& s : K )
    {
      if( s.dimension() != 1 )
        continue;

      auto u = s[0];         // vertex
      auto v = s[1];         // vertex

      writeEdge(out, u, v);
    }

    out << "\\end{tikzpicture}\n";
  }

  bool showVertexLabels() const noexcept
  {
    return _showVertexLabels;
  }

  void showVertexLabels( bool value ) noexcept
  {
    _showVertexLabels = value;
  }

  bool showBalls() const noexcept
  {
    return _showBalls;
  }

  void showBalls( bool value ) noexcept
  {
    _showBalls = value;
  }

  double ballRadius() const noexcept
  {
    return _ballRadius;
  }

  void ballRadius( double radius )
  {
    _ballRadius = radius;
  }

private:

  /**
    Auxiliary function for creating a point in TikZ format. The point
    will be shown as a circle of a given size.

    @param out Output stream to which the point will be appended
    @param v   Vertex index
    @param p   Coordinates for the vertex (only two dimensions will be used)
  */

  template <class Coordinate, class Index> void writePoint( std::ostream& out, Index v, Coordinate& p )
  {
    auto x = p[0];
    auto y = p[1];

    out << "\\coordinate";

    if( _showVertexLabels )
      out << "[label=" << _labelPosition << ":" << std::to_string(v) << "] ";

    out << "(" << v << ") at (" << x << "," << y << ");\n";

    out << "\\filldraw[" << _pointColour << "]"
        << " " << "(" << v << ") circle (" << _pointSize << _pointSizeUnit << ");\n";

    if( _showBalls )
    {
      out << "\\fill[" << _ballColour << ","
          << " fill opacity=" << _ballOpacity
          << "]"
          << " " << "(" << v << ") circle (" << _ballRadius << "cm" << ");\n";
    }
  }

  /**
    Auxiliary function for creating an edge in TikZ format. The edge
    will be shown as a line connecting the two vertices.

    @param out Output stream to which the edge will be appended
    @param u   Source vertex index
    @param v   Target vertex index
  */

  template <class Index> void writeEdge( std::ostream& out, Index u, Index v )
  {
    out << "\\draw[" << _lineColour << ","
        << " line width=" << _lineWidth << _lineWidthUnit
        <<"]"
        << " " << "(" << u << ") -- (" << v << ");\n";
  }

  bool _showBalls         = false;
  double _ballOpacity     = 0.1;
  double _ballRadius      = 0.0;
  std::string _ballColour = "black";

  bool _showVertexLabels     = false;
  std::string _labelPosition = "above";

  // Node/vertex configuration options ---------------------------------

  std::string _pointColour   = "black";
  std::string _pointSizeUnit = "pt";
  double _pointSize          = 1;

  // Line/edge configuration options -----------------------------------

  std::string _lineColour    = "black";
  std::string _lineWidthUnit = "mm";
  double _lineWidth          = 0.50;
};

} // namespace io

} // namespace topology

} // namespace aleph


#endif
