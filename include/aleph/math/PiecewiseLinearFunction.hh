#ifndef ALEPH_MATH_PIECEWISE_LINEAR_FUNCTION_HH__
#define ALEPH_MATH_PIECEWISE_LINEAR_FUNCTION_HH__

#include <aleph/math/KahanSummation.hh>

#include <functional>
#include <iterator>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <vector>

#include <cassert>
#include <cmath>

namespace aleph
{

namespace math
{

namespace detail
{

/** Performs linear interpolation between two points */
template <class D, class I> I lerp( D x, D x0, I y0, D x1, I y1 )
{
  return y0 + (y1-y0) * (x-x0) / (x1-x0);
}

} // namespace detail

/**
  @class PiecewiseLinearFunction
  @brief Models a piecewise linear function

  A piecewise linear function is fully defined by a set of pairs of
  coordinates, specifying values in the *domain* and the *image* of
  the function. This class permits various arithmetical operations,
  such as addition and subtraction, with piecewise linear functions
  and provides a number of other common operations.

  @tparam D Type of the *domain* of the function
  @tparam I Type of the *image* of the function
*/

template <class D, class I = D> class PiecewiseLinearFunction
{
public:
  using Domain = D;
  using Image  = I;

  /** Creates an empty piecewise linear function */
  PiecewiseLinearFunction() = default;

  /**
    Creates a new piecewise linear function from a range of values. The values
    must consist of pairs that are convertible to double. Duplicates are not
    permitted and will result in an exception being thrown.

    @param begin Input iterator to begin of range
    @param end   Input iterator to end of range
  */

  template <class InputIterator> PiecewiseLinearFunction( InputIterator begin,
                                                          InputIterator end )
  {
    for( InputIterator it = begin; it != end; ++it )
    {
      auto x       = it->first;
      auto y       = it->second;
      bool success = false;

      std::tie( std::ignore, success )
          = _data.insert( std::make_pair(x,y) );

      if( !success )
        throw std::runtime_error( "Duplicate value pairs not permitted for piecewise linear functions" );
    }

    // Inserts new points whenever the function intersects the x-axis in
    // order to ensure that the integration procedure works.
    if( !_data.empty() )
      this->insertIntersectionPoints();
  }

  // Evaluation --------------------------------------------------------

  /**
    Evaluates the piecewise linear function at a certain position. The
    function must not be evaluated beyond its domain. This will result
    in zeroes. For every other evaluation point, the function performs
    interpolation between the nearest values.

    @param x Coordinate at which to evaluate the function

    @returns Interpolated y value
  */

  Image operator()( Domain x ) const noexcept
  {
    auto range = _data.equal_range( x );
    auto begin = range.first;
    auto end   = range.second;

    // Beyond the domain
    if( begin == _data.end() || end == _data.begin() )
      return Image();

    // One of the stored values
    else if( begin->first == x )
      return begin->second;

    // Interpolation required
    auto right = begin;
    auto left  = std::prev( begin );

    return detail::lerp( x, left->first, left->second, right->first, right->second );
  }

  // Operations --------------------------------------------------------

  /** Calculates the sum of two piecewise linear functions */
  PiecewiseLinearFunction& operator+=( const PiecewiseLinearFunction& rhs ) noexcept
  {
    return this->apply( rhs, std::plus<Image>() );
  }

  /** Calculates the sum of two piecewise linear functions */
  PiecewiseLinearFunction operator+( const PiecewiseLinearFunction& rhs ) const noexcept
  {
    auto lhs = *this;
    lhs += rhs;
    return lhs;
  }

  /** Adds a given value to all points in the image of the function */
  PiecewiseLinearFunction operator+( Image lambda ) const noexcept
  {
    auto f = *this;
    return f.apply( lambda, std::plus<Image>() );
  }

  /** Adds a given value to all points in the image of the function */
  PiecewiseLinearFunction& operator+=( Image lambda ) const noexcept
  {
    return this->apply( lambda, std::plus<Image>() );
  }

  /** Calculates the difference of two piecewise linear functions */
  PiecewiseLinearFunction& operator-=( const PiecewiseLinearFunction& rhs ) noexcept
  {
    return this->apply( rhs, std::minus<Image>() );
  }

  /** Calculates the difference of two piecewise linear functions */
  PiecewiseLinearFunction operator-( const PiecewiseLinearFunction& rhs ) const noexcept
  {
    auto lhs = *this;
    lhs -= rhs;
    return lhs;
  }

  /** Unary minus: negates all values in the image of the piecewise linear function */
  PiecewiseLinearFunction operator-() const noexcept
  {
    PiecewiseLinearFunction f;

    for( auto&& pair : _data )
      f._data.insert( std::make_pair( pair.first, -pair.second ) );

    return f;
  }

  /** Subtracts a given value from all points in the image of the function */
  PiecewiseLinearFunction operator-( Image lambda ) const noexcept
  {
    auto f = *this;
    return f.apply( lambda, std::minus<Image>() );
  }

  /** Subtracts a given value from all points in the image of the function */
  PiecewiseLinearFunction& operator-=( Image lambda ) const noexcept
  {
    return this->apply( lambda, std::minus<Image>() );
  }

  /** Multiplies the given piecewise linear function with a scalar value */
  PiecewiseLinearFunction& operator*=( Image lambda ) noexcept
  {
    for( auto&& pair : _data )
      pair.second *= lambda;

    return *this;
  }

  /** Multiplies the given piecewise linear function with a scalar value */
  PiecewiseLinearFunction operator*( Image lambda ) const noexcept
  {
    auto f = *this;
    f *= lambda;
    return f;
  }

  /** Divides the given step function by a scalar value */
  PiecewiseLinearFunction& operator/=( I lambda )
  {
    if( lambda == I() )
      throw std::runtime_error( "Attempted division by zero" );

    return this->operator*=( 1/lambda );
  }

  /** Divides the given step function by a scalar value */
  PiecewiseLinearFunction operator/( I lambda ) const noexcept
  {
    auto f = *this;
    f /= lambda;
    return f;
  }

  // Queries -----------------------------------------------------------

  /** Copies the domain values to an output iterator */
  template <class OutputIterator> void domain( OutputIterator result ) const noexcept
  {
    for( auto&& pair : _data )
      *result++ = pair.first;
  }

  /** Copies the image values to an output iterator */
  template <class OutputIterator> void image( OutputIterator result ) const noexcept
  {
    for( auto&& pair : _data )
      *result++ = pair.second;
  }

  /**
    Compares two piecewise linear functions. The functions are
    considered to be *equal* when the take the same values for
    the same interval. This method will evaluate the functions
    over *all* points of the domain.

    @param rhs Other function to compare against

    @returns true if both functions are equal, regardless of whether one
    of the functions has a subdivided range of values.
  */

  bool operator==( const PiecewiseLinearFunction& rhs ) const noexcept
  {
    std::set<Domain> domain;

    this->domain( std::inserter( domain, domain.begin() ) );
      rhs.domain( std::inserter( domain, domain.begin() ) );

    bool result = true;

    for( auto&& x : domain )
    {
      auto y0 = this->operator()( x );
      auto y1 = rhs( x );

      if( y0 != y1 )
        return false;
    }

    return result;
  }

  /** Negates the comparison between two functions */
  bool operator!=( const PiecewiseLinearFunction& rhs ) const noexcept
  {
    return !this->operator==( rhs );
  }

  // Transformations ---------------------------------------------------

  /** Calculates the absolute value of the function */
  PiecewiseLinearFunction& abs() noexcept
  {
    for( auto&& pair : _data )
      pair.second = std::abs( pair.second );

    return *this;
  }

  /** Calculates the maximum (supremum) of the function */
  Image max() const noexcept
  {
    if( _data.empty() )
      return Image();

    auto max = std::numeric_limits<Image>::lowest();

    for( auto&& pair : _data )
      max = std::max( max, pair.second );

    return max;
  }

  /** Calculates the supremum (maximum) of the function */
  Image sup() const noexcept
  {
    return this->max();
  }

  /**
    Calculates the integral over the (absolute value of the) function,
    raised to the \f$p\f$-th power.
  */

  Image integral( Image p = Image(1) ) const noexcept
  {
    if( _data.empty() )
      return Image();

    // The Kahan summation ensures that small parts of the integral do
    // not disappear amidst the remaining values.
    KahanSummation<Image> norm = Image();

    auto previous = _data.begin();
    auto current  = std::next( _data.begin() );

    while( current != _data.end() )
    {
      // Coefficients of the line that connect the previous point and the current
      // point.
      auto x1 = previous->first;
      auto y1 = previous->second;
      auto x2 = current->first;
      auto y2 = current->second;
      auto m  = (y2 - y1) / (x2 - x1);
      auto c  = y1 - m * x1;

      // Evaluator for the integral. This is an application of Cavalieri's
      // quadrature formula.
      auto evaluator = [&]( Image x )
      {
        // Horizontal segments
        if( m == Image() )
          return std::pow( c, p ) * x;

        // Regular lines
        else
          return std::pow( m*x + c, p+1 ) / ( m * (p+1) );
      };

      auto integral   = std::abs( evaluator( x2 ) - evaluator( x1 ) );
      norm           += integral;

      previous = current;
      ++current;
    }

    return std::pow( norm, 1/p );
  }

private:

  /**
    Applies a binary operation to the current piecewise linear function
    and a given scalar. This is tantamount to *broadcasting* the result
    to a function that is entirely constant over a range of values. The
    operation does not require conversions or merges of the *domain* or
    the *range* of the input because only a scalar is involved.
  */

  template <class BinaryOperation> PiecewiseLinearFunction& apply( Image lambda,
                                                                   BinaryOperation operation )
  {
    for( auto&& pair : _data )
      pair.second = operation( pair.second, lambda );

    return *this;
  }

  /**
    Applies a binary operation to the current piecewise linear function
    and another function. To this, the domains of *both* functions will
    be merged, and the operation will be applied to  their values, with
    a suitable interpolation scheme in place.

    @param other     Second piecewise linear function
    @param operation Operation to apply to both functions

    @returns Reference to modified piecewise linear function. The
    original piecewise linear function ceases to exist.
  */

  template <class BinaryOperation> PiecewiseLinearFunction& apply( const PiecewiseLinearFunction& other,
                                                                   BinaryOperation operation )
  {
    std::set<Domain> xValues;

    for( auto&& pair : _data )
      xValues.insert( pair.first );

    for( auto&& pair : other._data )
      xValues.insert( pair.first );

    // Oherwise, the loop below will turn in an infinite loop since or
    // the validity of the `current` iterator would have to be checked
    // before-hand.
    if( xValues.empty() )
      return *this;

    // Intersection handling. This is required to ensure that the combination of
    // the two functions contains shared segments.
    {
      // This closure checks whether two line segments intersect. If this is the
      // case, the intersection point needs to be stored as an additional point
      // in the set of values.
      auto intersection = []( Domain x0, Image y0, Domain x1, Image y1,
                              Domain x2, Image y2, Domain x3, Image y3 )
      {
        auto s1x = x1 - x0;
        auto s1y = y1 - y0;
        auto s2x = x3 - x2;
        auto s2y = y3 - y2;

        auto s = ( -s1y * (x0-x2) + s1x * (y0-y2) ) / ( -s2x * s1y + s1x * s2y );
        auto t = (  s2x * (y0-y2) - s2y * (x0-x2) ) / ( -s2x * s1y + s1x * s2y );

        if( s >= 0 && s <= 1 && t >= 0 && t <= 1 )
          return x0 + t * s1x;

        return x0;
      };

      std::set<Domain> intersections;

      auto current = xValues.begin();
      auto next    = std::next( current );

      for( ; next != xValues.end(); )
      {
        auto x0 = *current;
        auto x1 = *next;

        auto y0 = this->operator()( x0 );
        auto y1 = this->operator()( x1 );
        auto y2 = other( x0 );
        auto y3 = other( x1 );

        auto x = intersection( x0, y0, x1, y1, x0, y2, x1, y3 );

        if( x != x0 )
          intersections.insert( x );

        current = next;
        next    = std::next( current );
      }

      xValues.insert( intersections.begin(), intersections.end() );
    }

    // Apply the operation to all points -------------------------------

    std::map<Domain, Image> data;

    for( auto&& x : xValues )
    {
      auto y1 = this->operator()( x );
      auto y2 =            other( x );

      data.insert( std::make_pair( x, operation( y1, y2 ) ) );
    }

    _data.swap( data );
    return *this;
  }

  /**
    Checks the segments of the piecewise linear function for intersections with
    the x-axis. In these cases, a new point needs to be inserted into the list
    of function values. Else, the integration procedure will not work because a
    segment might be considered to have an area of zero. The following example
    illustrates this problem:

    \verbatim
     /\
    /  \
    ---------
        \  /
         \/
    \endverbatim

    If the intersection in the middle is not being counted as a point on the PL
    function, the area of the corresponding segment will be zero. This is not
    the desired and expected behaviour.
  */

  void insertIntersectionPoints()
  {
    auto current = _data.begin();
    auto next    = std::next( current );

    std::set<Domain> intersections;

    for( ; next != _data.end(); )
    {
      auto x0 = current->first;
      auto y0 = current->second;

      auto x1 = next->first;
      auto y1 = next->second;

      // We do not need to check the other cases. If either one of the values is
      // zero, we already have an intersection.
      if( y0 * y1 < Image() )
      {
        auto m = ( y1 - y0 ) / ( x1 - x0 );
        auto c = y0 - m * x0;
        auto x = - c / m;

        intersections.insert( x );
      }

      current = next;
      next    = std::next( current );
    }

    for( auto&& x : intersections )
      _data.insert( std::make_pair( x, Image() ) );
  }

  /** Maps values in the domain to values in the image */
  std::map<Domain, Image> _data;
};

} // namespace math

} // namespace aleph

/** Multiplies a given piecewise linear function by a scalar value */
template <class D, class I, class T> aleph::math::PiecewiseLinearFunction<D, I> operator*( T lambda, const aleph::math::PiecewiseLinearFunction<D, I>& f )
{
  return f * lambda;
}

/** Divides a given piecewise linear function by a scalar value */
template <class D, class I, class T> aleph::math::PiecewiseLinearFunction<D, I> operator/( T lambda, const aleph::math::PiecewiseLinearFunction<D, I>& f )
{
  return f / lambda;
}

/**
  Output operator of a piecewise linear function. Permits printing the
  function values to an *output stream* in order to permit their usage
  by other functions and/or programs.
*/

template <class D, class I> std::ostream& operator<<( std::ostream& o, const aleph::math::PiecewiseLinearFunction<D, I>& f )
{
  std::vector<D> X;
  std::vector<I> Y;

  f.domain( std::back_inserter(X) );
  f.image ( std::back_inserter(Y) );

  assert( X.size() == Y.size() );

  for( std::size_t i = 0; i < X.size() && i < Y.size(); i++ )
    o << X[i] << "\t" << Y[i] << "\n";

  return o;
}

#endif
