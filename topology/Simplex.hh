#ifndef ALEPH_TOPOLOGY_SIMPLEX_HH__
#define ALEPH_TOPOLOGY_SIMPLEX_HH__

#include <boost/functional/hash.hpp>

#include <boost/iterator/iterator_adaptor.hpp>
#include <boost/iterator/filter_iterator.hpp>

#include <algorithm>
#include <initializer_list>
#include <iosfwd>
#include <stdexcept>
#include <vector>

namespace aleph
{

namespace topology
{

template <
  class D,
  class V = unsigned short
>
class Simplex
{
public:

  // Aliases & declarations -------------------------------------------
  //
  // Note that these aliases follow the STL conventions in order to make it
  // easier to use the class with STL algorithms.

  using DataType                      = D;
  using VertexType                    = V;

  using data_type                     = DataType;
  using vertex_type                   = VertexType;

  using vertex_container_type         = std::vector<vertex_type>;
  using vertex_iterator               = typename vertex_container_type::iterator;
  using const_vertex_iterator         = typename vertex_container_type::const_iterator;
  using reverse_vertex_iterator       = typename vertex_container_type::reverse_iterator;
  using const_reverse_vertex_iterator = typename vertex_container_type::const_reverse_iterator;

  // I cannot describe the class inline because boost::iterator_adaptor expects
  // a _complete_ class.
  class boundary_iterator;

  // Constructors ------------------------------------------------------

  /** Creates an empty simplex */
  Simplex()
    : _data( DataType() )
  {
  }

  /**
    Creates a new 0-simplex from the given vertex.

    @param u    Vertex
    @param data Data to assign to simplex
  */

  Simplex( VertexType u, DataType data = DataType() )
    : _vertices( 1, u )
    , _data( data )
  {
  }

  /**
    Creates a new simplex from another simplex while setting the data for the
    new simplex. This is handy for copying the vertices but not the data of a
    given simplex, which occurs in cases where the data for the new simplex is
    \i not determined by the simplex that is to be copied.

    @param simplex Simplex to copy vertices from
    @param data    Data to assign new simplex
  */

  explicit Simplex( const Simplex<DataType, VertexType>& simplex, DataType data )
    : _vertices( simplex._vertices )
    , _data( data )
  {
  }

  /**
    Creates a new simplex from a range of vertices. The input iterators are
    supposed to belong to a range of vertices for the simplex. This range need
    not be ordered. Note that the simplex will \b not check for duplicate
    values.

    @param begin Iterator to begin of vertex range
    @param end   Iterator to end of vertex range
    @param data  Data to assign to simplex
  */

  template <class InputIterator>
  Simplex( InputIterator begin, InputIterator end,
           DataType data = DataType() )
    : _vertices( begin, end )
    , _data( data )
  {
    std::sort( _vertices.begin(), _vertices.end(), std::greater<VertexType>() );

    // Ensures that the simplex does not contain the same vertex
    // multiple times. This may result in the empty simplex, but
    // there is no way around that.
    _vertices.erase(
      std::unique( _vertices.begin(), _vertices.end() ),
      _vertices.end()
    );
  }

  /**
    Creates a new simplex from a range of vertices. The vertices are not
    assumed to be ordered. It must be possible for me to convert them to
    the vertex type of the simplex.

    @param vertices Vertices
    @param data     Data to assign to simplex
  */

  template <class Vertex>
  Simplex( const std::initializer_list<Vertex>& vertices,
           DataType data = DataType() )
    : Simplex( vertices.begin(), vertices.end(), data )
  {
  }

  // vertices ----------------------------------------------------------

  /** @returns Iterator to begin of simplex vertex range */
  const_vertex_iterator begin() const
  {
    return _vertices.begin();
  }

  /** @returns Iterator to end of simplex vertex range */
  const_vertex_iterator end() const
  {
    return _vertices.end();
  }

  /** @returns Reverse begin iterator of simplex vertex range */
  const_reverse_vertex_iterator rbegin() const
  {
    return _vertices.rbegin();
  }

  /** @returns Reverse end iterator of simplex vertex range */
  const_reverse_vertex_iterator rend() const
  {
    return _vertices.rend();
  }

  /**
    Checks whether the current simplex contains a given vertex. This is
    required for intersection queries, for example.

    @param vertex Vertex to search for

    @returns true if the simplex contains the current vertex at least once,
    else false.
  */

  bool contains( VertexType vertex ) const
  {
    return std::find( this->begin(), this->end(), vertex ) != this->end();
  }

  // boundary ----------------------------------------------------------

  /** @returns Boundary iterator to begin of boundary */
  boundary_iterator begin_boundary() const
  {
    // Check dimension directly in order to handle the empty simplex
    if( _vertices.empty() || _vertices.size() <= 1 )
      return this->end_boundary();

    return boundary_iterator( _vertices.begin(), _vertices);
  }

  /** @returns Boundary iterator to end of boundary */
  boundary_iterator end_boundary() const
  {
    return boundary_iterator( _vertices.end(), _vertices);
  }

  // Data --------------------------------------------------------------

  /**
    Assigns the simplex a new value for its data object. This function does not
    perform any sanity checks. The value is simply copied and stored.

    @param data Data to assign
  */

  void setData( DataType data = DataType() )
  {
    _data = data;
  }

  /** @returns Current value of simplex data object */
  DataType data() const
  {
    return _data;
  }

  // Attribute access --------------------------------------------------

  /** @returns true if the simplex is empty, i.e. it has no vertices */
  bool empty() const
  {
    return _vertices.empty();
  }

  /**
    @returns true if the simplex is valid, i.e. it is not empty. This function
    allows the simplex class to be used in expressions such as this one:

    @code
    Simplex<double> mySimplex( 2, 2.0 ); // 0-simplex
    if( mySimplex )
    {
      // Do stuff with a valid simplex.
    }
    @endcode
  */

  explicit operator bool() const
  {
    // A bit verbose but more readable :)
    return this->empty() ? false : true;
  }

  /**
    @returns Dimension of simplex

    @throws std::runtime_error if the dimension of the empty simplex is
    queried. If you do this, it's your own fault.
  */

  std::size_t dimension() const
  {
    if( _vertices.empty() )
      throw std::runtime_error( "Querying dimension of empty simplex" );
    else
      return _vertices.size() - 1;
  }

  /** @returns Number of vertices of the simplex */
  std::size_t size() const
  {
    return _vertices.size();
  }

  /**
    Returns a vertex (specified by an index) of the current simplex.

    @param   index Index of vertex in simplex
    @returns Vertex of simplex, specified by an index.
    @throws  std::out_of_range if the index is out of range.
  */

  VertexType operator[]( std::size_t index ) const
  {
    return _vertices.at( index );
  }

  // Comparison operators ----------------------------------------------

  /**
    Checks whether two simplices are equal. Two simplices are considered equal
    if their vertices are being equal. Simplex data is \b not checked by this
    function as this would complicate simplex queries in a simplicial complex.

    @param other Simplex to check for equality
    @returns true if the two simplices are equal (see above), else false.
  */

  bool operator==( const Simplex& other ) const
  {
    return this->_vertices == other._vertices;
  }

  /**
    Checks whether two simplices are inequal. This function is simply the
    negation of operator==.

    @param other Simplex to check for inequality
    @returns true if the two simplices differ, else false.
  */

  bool operator!=( const Simplex& other ) const
  {
    return !this->operator==( other );
  }

  /**
    Comparison operator for simplices. Uses lexicographical comparison to
    obtain a weak ordering of the simplices. This is used in many algorithms.

    @param other Simplex to compare current simplex to
    @returns  true if current simplex is to be sorted before other simplex.
  */

  bool operator<( const Simplex& other ) const
  {
    return std::lexicographical_compare( this->_vertices.begin(), this->_vertices.end(),
                                         other._vertices.begin(), other._vertices.end() );
  }

  // Convenience functions ---------------------------------------------

  template <class DataType>                   friend std::size_t hash_value( const Simplex<DataType>& s );
  template <class DataType, class VertexType> friend std::size_t hash_value( const Simplex<DataType, VertexType>& s );

private:

  /**
    The vertices making up the current simplex. Their values are irrelevant,
    though it is assumed that the vertices represent some data points stored
    _outside_ the simplex.
  */

  vertex_container_type _vertices;

  /**
    Data stored within the simplex. The type and semantics of this member
    variable depend on the template parameters of the simplex class. For
    example, a \c double value is used at several places in order to add a \i
    weight or a \i distance for the simplex.
  */

  DataType _data;
};

// ---------------------------------------------------------------------

template <class DataType, class VertexType>
std::size_t hash_value( const Simplex<DataType, VertexType>& s )
{
  boost::hash< typename Simplex<DataType, VertexType>::vertex_container_type > hasher;
  return hasher( s._vertices );
}

template <class DataType> std::size_t hash_value( const Simplex<DataType>& s )
{
  return hash_value<DataType, typename Simplex<DataType>::vertex_type>( s );
}

// ---------------------------------------------------------------------

/**
  @class boundary_iterator
  @brief Iterator for traversing the boundary of a given simplex

  This iterator, inspired by Dmitriy Morozov's "Dionysus" framework, calculates
  the boundary of a given simplex while traversing it. Since the boundary
  simplices are created from scratch, they will _not_ have the correct weights
  set. Hence, the boundary iterator is only useful if some kind of lookup of
  generated simplices exists---as is the case for a simplicial complex, for
  example.
*/

template <
    class DataType,
    class VertexType
>
class Simplex<DataType, VertexType>::boundary_iterator
  : public boost::iterator_adaptor<boundary_iterator,
                                   const_vertex_iterator,
                                   Simplex<DataType, VertexType>,
                                   boost::use_default,
                                   Simplex<DataType, VertexType> >
{
public:

  using Iterator = const_vertex_iterator ;
  using Parent   = boost::iterator_adaptor<boundary_iterator,
                                           Iterator,
                                           Simplex<DataType, VertexType>,
                                           boost::use_default,
                                           Simplex<DataType, VertexType> >;

  /**
    Creates a new boundary iterator from a parent iterator (i.e. a simplex) and a
    set of vertices (i.e. the vertices of the parent simplex).

    @param it       Parent iterator
    @param vertices VertexType set
  */

  explicit boundary_iterator( Iterator it, const vertex_container_type& vertices )
    : Parent(it)
    , _vertices(vertices)
  {
  }

private:

  friend class boost::iterator_core_access;

  /** @returns Current boundary simplex */
  Simplex<DataType, VertexType> dereference() const
  {
    // This returns a new simplex. The simplex is created from a set of
    // vertices, which in turn is created by applying a filter to the set of
    // vertices stored in this iterator: Namely, the filter iterator will
    // return all vertices that are _not_ equal to its current position.

    vertex_container_type vertices(
          boost::make_filter_iterator( std::bind2nd( std::not_equal_to<vertex_type>(), *( this->base() ) ),
                                                     _vertices.begin(),
                                                     _vertices.end() ),
          boost::make_filter_iterator( std::bind2nd( std::not_equal_to<vertex_type>(), *( this->base() ) ),
                                                     _vertices.end(),
                                                     _vertices.end() )
          );

    return Simplex<DataType, VertexType>( vertices.begin(), vertices.end() );
  }

  /**
    Reference to vertex set; this is required because the boundary iterator
    iterates over a set of vertices and may thus not exist without one.
  */

  const vertex_container_type& _vertices;
};

// ---------------------------------------------------------------------

/**
  Outputs a simplex to an ostream. This is used for debugging purposes.

  @param o Output stream
  @param s Simplex to be added to o

  @returns Output stream with information about simplex s.
*/

template <class DataType, class VertexType>
std::ostream& operator<<( std::ostream& o, const topology::Simplex<DataType, VertexType>& s )
{
  auto numVertices = s.size();

  o << "{";

  for( decltype(numVertices) i = 0; i < numVertices; i++ )
  {
    if( i != 0 )
      o << " ";

    o << s[i];
  }

  if( s.data() != DataType() )
    o << " (" << s.data() << ")";

  o << "}";

  return o;
}

// ---------------------------------------------------------------------

}

}

#endif
