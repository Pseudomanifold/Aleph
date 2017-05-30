#ifndef ALEPH_TOPOLOGY_LOWER_STAR_HH__
#define ALEPH_TOPOLOGY_LOWER_STAR_HH__

#include <limits>
#include <vector>

namespace aleph
{

namespace filtrations
{

/**
  @class LowerStar
  @brief Lower-star filtration functor for simplicial complexes

  This functor calculates the lower-star filtration of a simplicial complex. To
  this end, it requires a vector containing the function values of each vertex
  in the simplicial complex. Since the functor will probably be \i copied a lot
  by the sorting procedure, it is advisable to wrap its usage in an std::ref
  object, like so:

  \code

    aleph::filtrations::LowerStar<T> functor( functionValues );

    simplicialComplex.sort( std::ref( functor ) );

  \endcode

  Otherwise, there might be spurious calls to the copy constructor, resulting
  in performance decreases.
*/

template <class Simplex> class LowerStar
{
public:

  using DataType   = typename Simplex::DataType;
  using VertexType = typename Simplex::VertexType;

  /**
    Creates a new lower-star filtration from a range of input iterators.

    @param begin Input iterator to begin of range
    @param end   Input iterator to end of range
  */

  template <class InputIterator> LowerStar( InputIterator begin,
                                            InputIterator end )
    : _values( begin, end )
  {
  }

  /**
    Using the function values stored by the functor, compares two simplices and
    checks which one comes first in a sorted sequence of simplices.

    @param s First simplex
    @param t Second simplex

    @returns true if simplex s is to be sorted before simplex t
  */

  bool operator()( const Simplex& s, const Simplex& t ) const
  {
    auto&& sValue = this->maximumValue( s );
    auto&& tValue = this->maximumValue( t );

    // Easy case: The values differ and s < t
    if( sValue < tValue )
      return true;

    // If the values of s and t coincide, the dimension decides about the order
    // of s and t. If their dimension also coincides, the tie is broken up by
    // using the lexicographical ordering.
    else if( sValue == tValue )
        return s < t;

    // Easy case: The values differ and s > t
    else
      return false;
  }

  /**
    Given a simplex, determines its maximum function value and returns the
    value.

    @param s Simplex
    @returns Maximum function value of the simplex
  */

  DataType maximumValue( const Simplex& s ) const
  {
    DataType maxValue = std::numeric_limits<DataType>::lowest();

    for( auto&& itVertex = s.begin(); itVertex != s.end(); ++itVertex )
      maxValue = std::max( maxValue, _values[ *itVertex ] );

    return maxValue;
  }

private:

  /** Stores function values for each vertex */
  std::vector<DataType> _values;
};

}

}

#endif
