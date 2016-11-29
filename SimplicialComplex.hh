#ifndef ALEPH_SIMPLICIAL_COMPLEX_HH__
#define ALEPH_SIMPLICIAL_COMPLEX_HH__

#include <boost/multi_index_container.hpp>

#include <boost/multi_index/indexed_by.hpp>
#include <boost/multi_index/mem_fun.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/random_access_index.hpp>

#include <algorithm>
#include <functional>
#include <initializer_list>
#include <iterator>
#include <iosfwd>
#include <limits>
#include <set>
#include <type_traits>
#include <vector>

namespace aleph
{

template <class Simplex> class SimplicialComplex
{
public:

  // Typedefs and tags for the multi-index container --------------------

  struct index_t {};
  struct lexicographical_t {};
  struct dimension_t {};

  using simplex_container_t
    =  boost::multi_index::multi_index_container<
        Simplex,
        boost::multi_index::indexed_by<
          boost::multi_index::random_access< boost::multi_index::tag<index_t> >,
          boost::multi_index::ordered_unique< boost::multi_index::tag<lexicographical_t>, boost::multi_index::identity<Simplex> >,
          boost::multi_index::ordered_non_unique< boost::multi_index::tag<dimension_t>, boost::multi_index::const_mem_fun<Simplex, std::size_t, &Simplex::dimension>        >
       >
  >;


  using const_iterator                 = typename simplex_container_t::template index<index_t>::type::const_iterator;
  using iterator                       = typename simplex_container_t::template index<index_t>::type::iterator;

  using lexicographical_iterator       = typename simplex_container_t::template index<lexicographical_t>::type::const_iterator;
  using const_lexicographical_iterator = typename simplex_container_t::template index<lexicographical_t>::type::iterator;

  using const_dimension_iterator       = typename simplex_container_t::template index<dimension_t>::type::const_iterator;
  using dimension_iterator             = typename simplex_container_t::template index<dimension_t>::type::iterator;

  // STL-like typedefs -------------------------------------------------

  // This simplifies writing algorithms that do not have any internal knowledge
  // about objects stored in this container.
  using value_type = Simplex;
  using ValueType  = value_type;

  // constructors ------------------------------------------------------

  /** Creates an empty simplicial complex. */
  SimplicialComplex()
  {
  }

  /**
    Creates a simplicial complex from an initializer list of simplices.

    @param simplices Simplices to insert into the simplicial complex
  */

  SimplicialComplex( std::initializer_list<Simplex> simplices )
    : _simplices( simplices.begin(), simplices.end() )
  {
    this->checkAndRestoreValidity();
  }

  /**
    Creates a simplicial complex from a given range of simplices. For
    flexibility reasons, the range type is left unspecified. Note that the
    complex will store _copies_ of the simplices for better memory management.

    @param begin  Iterator pointing to begin of range
    @param end    Iterator pointing to end of range
  */

  template <class InputIterator> SimplicialComplex( InputIterator begin, InputIterator end )
    : _simplices( begin, end )
  {
    this->checkAndRestoreValidity();
  }

  // simplex container modification ------------------------------------

  /** Clears the simplicial complex and removes all its simplices. */
  void clear()
  {
    _simplices.clear();
  }

  /**
    Given a range of simplices represented by two arbitrary input iterators,
    inserts the simplices into the simplicial complex.

    @param begin Iterator to begin of input range
    @param end   Iterator to end of input range
  */

  template <class InputIterator> void insert( InputIterator begin, InputIterator end )
  {
    _simplices.insert( _simplices.end(),
                       begin, end );

    this->checkAndRestoreValidity();
  }

  /**
    Given a range of simplices represented by two arbitrary input iterators,
    inserts the simplices into the simplicial complex \i without performing any
    validation.

    The same caveats as for the function push_back_without_validation() apply.

    @param begin Iterator to begin of input range
    @param end   Iterator to end of input range
  */

  template <class InputIterator> void insert_without_validation( InputIterator begin,
                                                                 InputIterator end )
  {
    _simplices.insert( _simplices.end(),
                       begin, end );
  }

  /**
    Inserts a new simplex into the simplicial complex. Note that the simplex is
    appended to the current filtration order, so the simplicial complex should
    be sorted again afterwards.

    @param simplex Simplex to insert into simplicial complex

    @warning If push_back() is called repeatedly, its performance will be worse
    than simply calling insert() for a range.
  */

  void push_back( const Simplex& simplex )
  {
    _simplices.push_back( simplex );

    this->checkAndRestoreValidity( simplex );
  }

  /**
    Inserts a new simplex into the simplicial complex \i without performing a
    validity check. By calling this function, a reduced simplicial complex that
    does not contain \i all subsets of every simplex can be created.

    Technically, the resulting simplicial complex ceases to be an abstract
    simplicial complex because it does not contain all possible faces. Note
    that the simplex is appended to the current filtration order, so the
    simplicial complex should be sorted again afterwards.

    @warning Other operations for modifying the simplicial complex might
    perform a validity check, thereby adding simplices (which you presumably
    want to avoid by calling this function). In order to avoid this problem,
    make sure to use this function \i exclusively for changing the simplicial
    complex.

    @param simplex Simplex to insert into simplicial complex
  */

  void push_back_without_validation( const Simplex& simplex )
  {
    _simplices.push_back( simplex );
  }

  /**
    Rearranges the simplices in the simplicial complex using an external view
    that contains each element exactly once. This function is useful when being
    applied to a view that has been modified by some STL algorithm (as the
    container only provides \i const iterators).

    @param first Input iterator to the beginning of the view
  */

  template <typename InputIterator> void rearrange( InputIterator first )
  {
    _simplices.rearrange( first );
  }

  /**
    Replaces a simplex stored in the simplicial complex (described by an
    iterator) by another simplex. Note that this function is the only way
    for modifying a fixed simplicial complex.

    @param position Iterator describing the simplex that is to be replaced
    @param simplex  Simplex to replace the simplex with

    @returns true if the replacement took place, else false.
  */

  bool replace( iterator position, const Simplex& simplex )
  {
    return _simplices.replace( position, simplex );
  }

  // Simplex container access ------------------------------------------

  /**
    @returns Iterator to begin of simplices in current filtration order. A
    filtration may be applied by sorting the simplicial complex.

    @see SimplicialComplex::sort()
  */

  iterator begin()
  {
    return _simplices.template get<index_t>().begin();
  }

  /** @overload begin() */
  const_iterator begin() const
  {
    return _simplices.template get<index_t>().begin();
  }

  /**
    @returns Iterator to end of simplices in current filtration order. A
    filtration may be applied by sorting the simplicial complex.

    @see SimplicialComplex::sort()
  */

  iterator end()
  {
    return _simplices.template get<index_t>().end();
  }

  /** @overload end() */
  const_iterator end() const
  {
    return _simplices.template get<index_t>().end();
  }

  /**
    @param   index Simplex index
    @returns Simplex at corresponding index position. Invalid indices will not
    be caught.
  */

  const Simplex& operator[]( std::size_t index ) const
  {
    return _simplices.template get<index_t>().operator[]( index );
  }

  /**
    @param   index Simplex index
    @returns Simplex at corresponding index position
    @throws  std::runtime_error for invalid indices
  */

  const Simplex& at( std::size_t index ) const
  {
    return _simplices.template get<index_t>().at( index );
  }

  /** @returns Iterator to begin of simplices in lexicographical order. */
  lexicographical_iterator begin_lexicographical()
  {
    return _simplices.template get<lexicographical_t>().begin();
  }

  /** @overload begin_lexicographical() */
  const_lexicographical_iterator begin_lexicographical() const
  {
    return _simplices.template get<lexicographical_t>().begin();
  }

  /** @returns Iterator to end of simplices in lexicographical order. */
  lexicographical_iterator end_lexicographical()
  {
    return _simplices.template get<lexicographical_t>().end();
  }

  /** @overload end_lexicographical() */
  const_lexicographical_iterator end_lexicographical() const
  {
    return _simplices.template get<lexicographical_t>().end();
  }

  /** @returns Iterator to begin of simplices in dimensional order. */
  dimension_iterator begin_dimension()
  {
    return _simplices.template get<dimension_t>().begin();
  }

  /** @overload begin_dimension() */
  const_dimension_iterator begin_dimension() const
  {
    return _simplices.template get<dimension_t>().begin();
  }

  /** @returns Iterator to end of simplices in dimensional order. */
  dimension_iterator end_dimension()
  {
    return _simplices.template get<dimension_t>().end();
  }

  /** @overload end_dimension() */
  const_dimension_iterator end_dimension() const
  {
    return _simplices.template get<dimension_t>().end();
  }

  /**
    Given an output iterator, calculates the vertex set of the simplicial
    complex. The vertex set contains all vertices that occur in at least
    one simplex that is stored in the simplicial complex. The vertices are
    guaranteed to be reported in ascending order.

    @param result Output iterator for storing the result
  */

  template <class OutputIterator> void vertices( OutputIterator result ) const
  {
    // Using a set as an intermediate storage has the advantage that all of the
    // vertices are guaranteed to be sorted, regardless of the current order of
    // the simplicial complex.

    std::set<typename Simplex::vertex_type> vertices;

    const_dimension_iterator d_it;
    const_dimension_iterator d_it_end;

    for( std::tie( d_it, d_it_end ) = this->range( 0 );
         d_it != d_it_end;
         ++d_it )
    {
      vertices.insert( *( d_it->begin() ) );
    }

    std::copy( vertices.begin(), vertices.end(), result );
  }

  /**
    Checks whether the simplicial complex contains a given simplex. Note that
    this only checks for vertex equality. The (optional) user data is not taken
    into account. This is by design, because a simplicial complex cannot
    contain two simplices with equal vertex sets, regardless of their user data
    values.

    @param simplex Simplex whose existence is checked
    @returns true if the simplicial complex contains the simplex, else false.
  */

  bool contains( const Simplex& simplex ) const
  {
    return _simplices.template get<lexicographical_t>().find( simplex ) != _simplices.template get<lexicographical_t>().end();
  }

  /**
    Searches the simplicial complex for a given simplex and, if found, returns
    an iterator to it. If the complex does not contain the simplex, an iterator
    to the end of the simplices is returned. No (optional) user data is taken
    into account here because this makes the query more flexible: The client
    may take \i any simplex as the input parameter --- without knowing the user
    data --- but the result simplex (if valid) contains known user data that
    may only be stored within the simplicial complex.

    @param simplex Simplex to query complex for

    @returns Iterator to simplex, or an iterator to the end of the simplicial
    complex if the complex does not contain the given simplex.
  */

  const_iterator find( const Simplex& simplex ) const
  {
    auto&& it
      = _simplices.template get<lexicographical_t>().find( simplex );

    if( it != this->end_lexicographical() )
      return _simplices.template project<index_t>( it );
    else
      return this->end();
  }

  /** @overload find() */
  iterator find( const Simplex& simplex )
  {
    auto&& it
      = _simplices.template get<lexicographical_t>().find( simplex );

    if( it != this->end_lexicographical() )
      return _simplices.template project<index_t>( it );
    else
      return this->end();
  }

  /**
    Given a simplex contained by the simplicial complex, looks up its index in
    the current filtration order. If the simplex is not part of the complex, an
    exception will be thrown.

    @param simplex Simplex

    @returns Index of simplex in current filtration

    @throws std::runtime_error if the simplex is not part of the simplicial
    complex.
  */

  std::size_t index( const Simplex& simplex ) const
  {
    auto&& itSimplex = this->find( simplex );

    if( itSimplex != this->end() )
      return static_cast<std::size_t>( std::distance( this->begin(), itSimplex ) );
    else
      throw std::runtime_error( "Queried simplex does not exist" );
  }

  /** @returns Number of simplices stored in simplicial complex */
  std::size_t size() const
  {
    return _simplices.size();
  }

  /**
    @returns true if the simplicial is empty, i.e. if it does not contain any
    simplices.
  */

  bool empty() const
  {
    return _simplices.empty();
  }

  /** @returns Maximum dimension of simplices stored in simplicial complex */
  std::size_t dimension() const
  {
    if( !this->empty() )
      return _simplices.template get<dimension_t>().rbegin()->dimension();
    else
      throw std::runtime_error( "Unable to query dimensionality of empty simplicial complex" );
  }

  // range queries -----------------------------------------------------

  /**
    Given a dimension, extracts all simplices whose dimension matches the
    user-specified one, and returns a pair of iterators for this range.

    @param dimension Dimension to extract simplices from

    @returns Pair of iterators describing the range of simplices matching the
    dimension. Note that the range is allowed to be empty.
  */

  std::pair<const_dimension_iterator, const_dimension_iterator> range( std::size_t dimension ) const
  {
    return this->range( [&] ( std::size_t d ) { return d >= dimension; },
                        [&] ( std::size_t d ) { return d <= dimension; } );
  }

  /**
    Given predicates describing the lower and upper bounds of a range, returns
    a pair of iterators for this range, assuming it is contained within the
    simplicial complex.

    @param lower Predicate describing lower bound of range
    @param upper Predicate describing upper bound of range

    @returns Pair of iterators describing the requested range.
   */

  template <typename LowerBounder, typename UpperBounder>
  std::pair<const_dimension_iterator, const_dimension_iterator> range( LowerBounder lower,
                                                                       UpperBounder upper ) const
  {
    return _simplices.template get<dimension_t>().range( lower, upper );
  }

  // filtration modification -------------------------------------------

  /**
    Allows changing the current order of simplices, i.e. applying a certain
    simplicial filtration. This function only requires a simplex comparison
    function as its input. See the Simplex class for admissable functors.

    @param comparison Simplex comparison object (or function)

    @see Simplex<T>::DataComparison
    @see Simplex<T>::DataDimensionComparison
  */

  template <class Comparison> void sort( Comparison&& comparison )
  {
    _simplices.sort( std::ref( comparison ) );
  }

  // -------------------------------------------------------------------

  /**
    Uses a range of vertex weights to recalculate all weights in the simplicial
    complex. Each higher-dimensional simplex is assigned the maximum of the
    weights of its lower-dimensional faces.

    @param begin Input iterator to begin of range
    @param end   Input iterator to end of range
  */

  template <class InputIterator> void recalculateWeights( InputIterator begin,
                                                          InputIterator end )
  {
    using data_type_  = typename std::iterator_traits<InputIterator>::value_type;
    using data_type   = typename Simplex::data_type;
    using vertex_type = typename Simplex::vertex_type;

    static_assert( std::is_same<data_type_, data_type>::value, "Data types must agree" );

    std::vector<data_type> weights( begin, end );

    // Reset all weights -----------------------------------------------

    for( auto itSimplex = this->begin(); itSimplex != this->end(); ++itSimplex )
    {
      Simplex newSimplex( *itSimplex );
      newSimplex.setData( std::numeric_limits<data_type>::max() );

      this->replace( itSimplex, newSimplex );
    }

    // Assign 0-dimensional weights ------------------------------------

    for( auto itSimplex = this->begin_dimension();
         itSimplex != this->end_dimension();
         ++itSimplex )
    {
      if( itSimplex->dimension() == 0 )
      {
        vertex_type v = *( itSimplex->begin() );

        Simplex newSimplex = ( *itSimplex );
        newSimplex.setData( weights.at(v) );

        this->replace( _simplices.template project<index_t>( itSimplex ),
                       newSimplex );
      }
      else
        break;
    }

    this->recalculateWeights();
  }

  // -------------------------------------------------------------------

  /**
    Recalculates simplex weights by assigning each simplex the maximum weight
    of its faces.

    @param skipOneDimensionalSimplices If set, skips both 0-dimensional and
    1-dimensional simplices and accepts their weights as the given truth.
  */

  void recalculateWeights( bool skipOneDimensionalSimplices = false )
  {
    // Assign weights based on lower-dimensional simplices -------------

    for( auto itSimplex = this->begin_dimension();
         itSimplex != this->end_dimension();
         ++itSimplex )
    {
      if(    ( itSimplex->dimension() == 0 )
          || ( skipOneDimensionalSimplices && itSimplex->dimension() == 1 ) )
      {
        continue;
      }

      typename Simplex::data_type weight
        = std::numeric_limits<typename Simplex::data_type>::lowest();

      for( auto itBoundary = itSimplex->begin_boundary();
           itBoundary != itSimplex->end_boundary();
           ++itBoundary )
      {
        const_iterator itPos = this->find( *itBoundary );
        if( itPos != this->end() )
        {
          weight = std::max( weight,
                             itPos->data() );
        }

        // The if-branch above ignores missing boundaries. This is useful when a
        // filtration is only partially defined, e.g. only up to the 2-simplices
        // and higher-dimensional simplices.
      }

      Simplex newSimplex = ( *itSimplex );
      newSimplex.setData( weight );

      this->replace( _simplices.template project<index_t>( itSimplex ),
                     newSimplex );
    }
  }

  // Container modification --------------------------------------------

  /**
    Allows simplex removal by value. The simplicial complex will check whether
    the given simplex exists. If so, it will be erased. Note that erasing a
    simplex may trigger the removal of co-faces of the simplex in order to
    remain validity.

    @param simplex Simplex to remove
  */

  void remove( const Simplex& simplex )
  {
    _simplices.template get<index_t>().remove( simplex );

    bool foundInvalidSimplex = false;

    do
    {
      foundInvalidSimplex = false;

      for( auto itSimplex = this->begin();
           itSimplex != this->end(); )
      {
        if( this->checkValidity( *itSimplex ) == false )
        {
          foundInvalidSimplex = true;
          itSimplex           = _simplices.template get<index_t>().erase( itSimplex );
        }
        else
          ++itSimplex;
      }
    }
    while( foundInvalidSimplex );
  }

  // comparison --------------------------------------------------------

  /**
    Checks two simplicial complexes for equality. This operator will
    check simplicial complexes for equality but only with respect to
    their current filtration order.
  */

  bool operator==( const SimplicialComplex<Simplex>& other ) const
  {
    if( this->size() != other.size() )
      return false;

    for( auto it1 = this->begin(), it2 = other.begin();
         it1 != this->end() && it2 != other.end();
         ++it1, ++it2 )
    {
      if( *it1 != *it2 )
        return false;
    }

    return true;
  }

  /**
    Checks whether two simplicial complexes differ by at least one
    simplex with respect to their current filtration order.
  */

  bool operator!=( const SimplicialComplex<Simplex>& other ) const
  {
    return this->operator==( other );
  }

  // debug -------------------------------------------------------------

  /**
    Adds information about a simplicial complex to an ostream. This function is
    useful for debug purposes or intensive logging.

    @param o ostream to add simplicial complex to
    @param S Simplicial complex to stream to ostream

    @returns ostream with information about simplicial complex
  */

  template <class Simplex_> friend std::ostream& operator<<( std::ostream& o,
                                                             const SimplicialComplex< Simplex_ >& S );

private:

  /**
    Checks and restores validity of the simplicial complex. A simplicial
    complex is deemed valid if for every simplex s, it contains every face f of
    s.
  */

  void checkAndRestoreValidity()
  {
    // If I encounter a face that is not stored in the simplicial complex, I
    // create the corresponding simplex and add it. This ensures that it is
    // possible to construct a simplicial complex from "partial" ranges of
    // simplices, e.g. a list of high-dimensional simplices whose faces need to
    // be calculated manually, for example.

    for( auto itSimplex = this->begin(); itSimplex != this->end(); itSimplex++ )
      this->checkAndRestoreValidity( *itSimplex );
  }

  /**
    Checks and restores validity of the simplicial complex after adding a
    single simplex. This means that the simplicial complex will check whether
    it contains \i all faces of the given simplex, as well as the simplex
    itself.

    If a face is found to be missing, it will be created by the function.

    @param simplex Simplex whose addition to the simplicial complex caused the
    validity check to be triggered.
  */

  void checkAndRestoreValidity( const Simplex& simplex )
  {
    for( auto itFace = simplex.begin_boundary();
         itFace != simplex.end_boundary();
         ++itFace )
    {
      const_lexicographical_iterator itStoredFace
          = _simplices.template get<lexicographical_t>().find( *itFace );

      if( itStoredFace == _simplices.template get<lexicographical_t>().end() )
      {
        // The new simplex shall contain the same vertices as the "face
        // simplex", but the data from its parent simplex. This ensures that
        // the data of a coface is always greater than or equal to the data of
        // its faces (assuming that the data type is comparable).

        _simplices.push_back( Simplex( *itFace,
                                       simplex.data() ) );
      }
    }
  }

  /**
    Checks validity of a single simplex. A simplex in the simplicial complex is
    deemed valid if all of its faces can be found in the complex.

    @param simplex Simplex whose validity is to be checked
    @returns true if the simplex is valid, else false
  */

  bool checkValidity( const Simplex& simplex )
  {
    for( auto itFace = simplex.begin_boundary();
         itFace != simplex.end_boundary();
         ++itFace )
    {
      // Check whether an unknown face has been found. If so, the simplex is
      // invalid.
      const_lexicographical_iterator itStoredFace
          = _simplices.template get<lexicographical_t>().find( *itFace );

      if( itStoredFace == _simplices.template get<lexicographical_t>().end() )
        return false;
    }

    return true;
  }

  /**
    Simplex container. boost::multi_index is used to provide different "views"
    to the simplicial data set. The first view uses the current sorting order
    of simplices. It corresponds to a given "filtration" order, in which faces
    precede cofaces. The second view allows a lexicographical traversal of the
    current simplicial complex. This view is used to provide fast index
    queries.
  */

  simplex_container_t _simplices;
};

// ---------------------------------------------------------------------

template <class Simplex> std::ostream& operator<<( std::ostream& o,
                                                   const SimplicialComplex<Simplex>& S )
{
  if( S.empty() )
    return o;

  o << std::string( 80, '-' ) << "\n";

  for( auto it = S.begin(); it != S.end(); ++it )
    o << *it << "\n";

  o << std::string( 80, '-' ) << "\n";

  return o;
}

// ---------------------------------------------------------------------

}

#endif
