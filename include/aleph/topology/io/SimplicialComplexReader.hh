#ifndef ALEPH_TOPOLOGY_IO_SIMPLICIAL_COMPLEX_READER_HH__
#define ALEPH_TOPOLOGY_IO_SIMPLICIAL_COMPLEX_READER_HH__

#include <algorithm>
#include <fstream>
#include <iterator>
#include <string>
#include <stdexcept>
#include <vector>

#include <aleph/topology/io/EdgeLists.hh>
#include <aleph/topology/io/GML.hh>
#include <aleph/topology/io/HDF5.hh>
#include <aleph/topology/io/Pajek.hh>
#include <aleph/topology/io/PLY.hh>
#include <aleph/topology/io/VTK.hh>

#include <aleph/utilities/Filesystem.hh>

namespace aleph
{

namespace topology
{

namespace io
{

/**
  @class SimplicialComplexReader
  @brief Generic simplicial complex reader class

  The purpose of this class is to provide a unified interface for
  reading a simplicial complex from an input file, assign weights
  that are consistent, and sort it.

  The class offers a rich interface but not every file format has
  the proper capabilities to support it.
*/

class SimplicialComplexReader
{
public:

  /**
    Attempts to read the simplicial complex \p K from the given input
    file, while using a default strategy for assigning weights. For a
    support file format, weights of higher-dimensional simplices will
    be assigned according to the *maximum* of their vertices.

    @param filename Input file
    @param K        Simplicial complex to read

    @see SimplicialComplexReader::operator()( const std::string&, SimplicialComplex&, Functor )
  */

  template <class SimplicialComplex> void operator()( const std::string& filename, SimplicialComplex& K )
  {
    using Simplex  = typename SimplicialComplex::ValueType;
    using DataType = typename Simplex::DataType;

    this->operator()( filename, K, [] ( DataType a, DataType b ) { return std::max(a,b); } );
  }

  /**
    Reads a simplicial complex from a file using a separate functor to
    assign weights for higher-dimensional simplices. The functor needs
    to support the following interface:

    \code{.cpp}
    using SimplexType = typename SimplicialComplex::ValueType;
    using DataType    = typename Simplex::DataType;

    DataType Functor::operator()( DataType a, DataType b )
    {
      // Do something with a and b. Typically, this would be calculating
      // either the minimum or the maximum...
      return std::max(a, b);
    }
    \endcode

    Typically, if you use `std::max()` in the functor, you obtain
    a simplicial complex that is filtrated by its *sublevel* sets,
    whereas you obtain a *superlevel* set filtration whenever you
    use `std::min()` instead. Other functors need to yield a valid
    filtration. This class does *not* check for it!

    @see SimplicialComplexReader::operator()( const std::string&, SimplicialComplex& )
  */

  template <class SimplicialComplex, class Functor> void operator()( const std::string& filename, SimplicialComplex& K, Functor functor )
  {
    std::ifstream in( filename );
    if( !in )
      throw std::runtime_error( "Unable to read input file" );

    auto extension = aleph::utilities::extension( filename );

    // The GML parser works more or less on its own and does not make
    // use of any stored variables.
    if( extension == ".gml" )
    {
      GMLReader reader;
      reader( filename, K );

      using VertexType = typename SimplicialComplex::ValueType::VertexType;

      _labels = getLabels( reader.getNodeAttribute( _labelAttribute ), reader.id_to_index<VertexType>() );
    }

    // The HDF5 parser permits the use of a functor that assigns
    // a weight to a higher-dimensional simplex.
    else if( extension == ".h5" )
    {
      HDF5SimpleDataSpaceReader reader;
      reader( filename, K );
    }

    // The Pajek parser also works on its own and does not require any
    // other configuration options.
    else if( extension == ".net" )
    {
      PajekReader reader;
      reader( filename, K );

      _labels = getLabels( reader.getLabelMap() );
    }

    // For PLY files, the parser supports optionally reading a property
    // for every vertex to assign as a data property to the vertices of
    // the simplicial complex.
    else if( extension == ".ply" )
    {
      PLYReader reader;

      if( !_dataAttribute.empty() )
        reader.setDataProperty( _dataAttribute );

      reader( filename, K );
    }

    // VTK files permit the use of a functor that assigns a weight to a
    // higher-dimensional simplex.
    else if( extension == ".vtk" )
    {
      VTKStructuredGridReader reader;
      reader( filename, K, functor );
    }

    // In all other cases, we fall back to reading a graph from an edge
    // list, with optional weights being specified.
    //
    // TODO: does it make sense to make this reader configurable?
    else
    {
      EdgeListReader reader;
      reader.setTrimLines();
      reader.setReadWeights();

      reader( filename, K );
    }
  }

  /**
    Sets the attribute that is used to extract data values from input
    files. For PLY files, for example, this means using an  attribute
    of the vertices of the mesh. Note that the attribute is only used
    if non-empty.
  */

  void setDataAttribute( const std::string& attribute ) noexcept
  {
    _dataAttribute = attribute;
  }

  /**
    @returns Current data attribute
    @see SimplicialComplexReader::setDataAttribute()
  */

  std::string dataAttribute() const noexcept
  {
    return _dataAttribute;
  }

  /** Sets current label attribute */
  void setLabelAttribute( const std::string& attribute ) noexcept
  {
    _labelAttribute = attribute;
  }

  /** @returns Current label attribute */
  std::string labelAttribute() const noexcept
  {
    return _labelAttribute;
  }

  /**
    @returns Node labels (if any). The indexing follows the vertex order
    in the simplicial complex and the graph. The label at index `0` thus
    corresponds to the vertex `0` in the resulting simplicial complex.
  */

  std::vector<std::string> labels() const noexcept
  {
    return _labels;
  }

private:

  /**
    Specifies an attribute of the input data file that stores the data
    property of the corresponding simplicial complex. The usage of the
    attribute depends on the file format. For a PLY file, for example,
    this means that a certain vertex attribute of the mesh set is used
    to obtain weights.

    This property is only used if it is not left empty.

    @see SimplicialComplexReader::dataAttribute()
    @see SimplicialComplexReader::setDataAttribute()
  */

  std::string _dataAttribute;

  /**
    Specifies an attribute of the input data file that stores labels for
    the corresponding simplicial complex. Since not all file formats are
    capable of storing labels at all, this attribute may not be used.

    @see SimplicialComplexReader::labelAttribute()
    @see SimplicialComplexReader::setLabelAttribute()
  */

  std::string _labelAttribute = "label";

  /**
    Converts a map of labels, obtained from a subordinate reader class,
    to a vector of labels. The order is guaranteed to follow the vertex
    order in the graph and in the simplicial complex.
  */

  template <class Map> std::vector<std::string> getLabels( const Map& map )
  {
    if( map.empty() )
      return {};

    std::vector<std::string> labels;
    labels.reserve( map.size() );

    for( auto&& pair : map )
    {
      if( !pair.second.empty() )
        labels.push_back( pair.second );
    }

    return labels;
  }

  // TODO: missing documentation
  template <class Map1, class Map2> std::vector<std::string> getLabels( const Map1& labelMap, const Map2& idMap )
  {
    if( labelMap.empty() || idMap.empty() )
      return {};

    std::vector<std::string> labels;
    labels.resize( labelMap.size() );

    for( auto&& pair : labelMap )
    {
      if( !pair.second.empty() )
        labels.at( idMap.at( pair.first ) ) =  pair.second;
    }

    return labels;
  }

  /**
    Optionally stores labels that have been extracted when reading an
    input file.

    Use SimplicialComplexReader::labels() to access them.
  */

  std::vector<std::string> _labels;
};

} // namespace io

} // namespace topology

} // namespace aleph

#endif
