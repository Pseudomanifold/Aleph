#ifndef ALEPH_UPPER_STAR_HH__
#define ALEPH_UPPER_STAR_HH__

#include <limits>
#include <vector>

namespace aleph
{

namespace filtrations
{

/**
  @class UpperStar
  @brief Upper-star filtration functor for simplicial complexes

  This functor calculates the upper-star filtration of a simplicial complex. To
  this end, it requires a vector containing the function values of each vertex
  in the simplicial complex.

  @see LowerStar
*/

template <class Simplex> class UpperStar
{
public:

  using DataType   = typename Simplex::DataType;
  using VertexType = typename Simplex::VertexType;

  /**
    Creates a new upper-star filtration from a range of input iterators.

    @param begin Input iterator to begin of range
    @param end   Input iterator to end of range
  */

  template <class InputIterator> UpperStar( InputIterator begin, InputIterator end )
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
    auto&& sValue = this->minimumValue( s );
    auto&& tValue = this->minimumValue( t );

    // Easy case: The values differ and s > t
    if( sValue > tValue )
      return true;

    // If the values of s and t coincide, the dimension decides about their
    // order. If their dimension also coincides, the tie is broken up by using
    // the lexicographical ordering.
    else if( sValue == tValue )
      return s < t;

    // Easy case: The values differ and s > t
    else
      return false;
  }

  /**
    Given a simplex, determines its minimum function value and returns the
    value.

    @param s Simplex
    @returns Minimum function value of the simplex
  */

  DataType minimumValue( const Simplex& s ) const
  {
    DataType minValue = std::numeric_limits<DataType>::max();

    for( auto&& itVertex = s.begin(); itVertex != s.end(); ++itVertex )
      minValue = std::min( minValue, _values[ *itVertex ] );

    return minValue;
  }

private:

  /** Stores function values for each vertex */
  std::vector<DataType> _values;
};

}

}

#endif
