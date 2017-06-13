#ifndef ALEPH_MATH_STEP_FUNCTION_HH__
#define ALEPH_MATH_STEP_FUNCTION_HH__

#include <iterator>
#include <limits>
#include <ostream>
#include <set>
#include <stdexcept>
#include <vector>

#include <cmath>
#include <cstdlib>

namespace aleph
{

namespace math
{

namespace detail
{

template <class T> T next( T x )
{
  if( std::numeric_limits<T>::is_integer )
    return x+1;
  else
    return std::nextafter( x, std::numeric_limits<T>::max() );
}

template <class T> T previous( T x )
{
  if( std::numeric_limits<T>::is_integer )
    return x-1;
  else
    return std::nextafter( x, std::numeric_limits<T>::lowest() );
}

} // namespace detail

template <class D, class I = D> class StepFunction
{
public:

  using Domain = D;
  using Image  = I;

  /**
    Auxiliary class for representing an indicator function interval of the step
    function. Each indicator function is only non-zero within its interval, and
    zero outside.
  */

  class IndicatorFunction
  {
  public:
    IndicatorFunction( D a )
      : _a( a    )
      , _b( a    )
      , _y( I(1) )
    {
    }

    IndicatorFunction( D a, D b, I y )
      : _a( a )
      , _b( b )
      , _y( y )
    {
      if( a > b )
        throw std::runtime_error( "Invalid interval specified" );
    }

    const D& a() const noexcept { return _a; }
    const D& b() const noexcept { return _b; }

          I& y()       noexcept { return _y; }
    const I& y() const noexcept { return _y; }

    bool contains( D x ) const noexcept
    {
      return this->a() <= x && x <= this->b();
    }

    /** Standard (signed) integral */
    I integral() const noexcept
    {
      return this->y() * static_cast<I>( ( this->b() - this->a() ) );
    }

    /** Unsigned integral raised to a certain power */
    I integral_p( I p ) const noexcept
    {
      auto value  = std::abs( this->integral() );
      return std::pow( value, p );
    }

    I operator()( D x ) const noexcept
    {
      if( this->contains( x ) )
        return this->y();
      else
        return I();
    }

    bool operator<( const IndicatorFunction& other ) const
    {
      // Permits that intervals intersect in a single point, as this is
      // simplifies the composition of multiple step functions.
      return this->b() <= other.a();
    }

  private:
    D _a;
    D _b;
    I _y;
  };

  /** Adds a new indicator function to the step function */
  void add( D a, D b, I y ) noexcept
  {
    _indicatorFunctions.insert( IndicatorFunction(a,b,y) );
  }

  /** Returns the domain of the function */
  template <class OutputIterator> void domain( OutputIterator result ) const noexcept
  {
    for( auto&& f : _indicatorFunctions )
    {
      *result++ = f.a();
      *result++ = f.b();
    }
  }

  /** Returns the image of the function */
  template <class OutputIterator> void image( OutputIterator result ) const noexcept
  {
    for( auto&& f : _indicatorFunctions )
      *result++ = f.y();
  }

  /** Returns the function value at a certain position */
  I operator()( D x ) const noexcept
  {
    I value = I();

    for( auto&& f : _indicatorFunctions )
    {
      // TODO: Not sure whether I really want this. The step functions must not
      // overlap anyway...
      if( f.contains(x) && std::abs( f(x) ) > value )
        value = f(x);
    }

    return value;
  }

  /** Calculates the sum of this step function with another step function */
  StepFunction operator+( const StepFunction& other ) const noexcept
  {
    auto&& f = *this;
    auto&& g = other;

    std::set<D> domain;

    f.domain( std::inserter( domain, domain.begin() ) );
    g.domain( std::inserter( domain, domain.begin() ) );

    StepFunction<D,I> h;

    if( domain.empty() )
      return h;

    auto prev = *domain.begin();
    auto curr = *domain.begin();

    for( auto it = std::next( domain.begin() ); it != domain.end(); ++it )
    {
      curr    = *it;

      auto y1 = f( prev );
      auto y2 = g( prev );
      auto y3 = f( curr );
      auto y4 = g( curr );
      auto y5 = f( detail::next( prev ) );
      auto y6 = g( detail::next( prev ) );

      if( y1 == y3 && y2 == y4 )
        h.add( prev, curr, y1+y2 );

      // A subdivision of the interval [prev, curr] is required because
      // at least one of the functions differs on the end points of the
      // interval.
      else if( y1 != y3 || y2 != y4 )
      {
        // Make the interval smaller if we hit a point where at least a
        // single function changes.
        if( y1 != y5 || y2 != y6 )
          prev = detail::next( prev );

        if( y5+y6 != y3+y4 )
        {
          // [prev, curr - epsilon]: y5+y6
          // [curr, curr + epsilon]: y3+y4
          h.add( prev, detail::previous( curr ), y5+y6 );
          h.add( curr, detail::next(     curr ), y3+y4 );
        }

        // FIXME: Not sure whether this is really the best workaround
        // here or whether I need a different strategy.
        else
          h.add( prev, detail::next( curr ), y3+y4 );

        // Ensures that the next interval uses the proper start point for the
        // indicator function interval.
        curr = detail::next( curr );
      }

      prev = curr;
    }

    return h;
  }

  /** Unary minus: negates all values in the image of the step function */
  StepFunction operator-() const noexcept
  {
    StepFunction f;

    for( auto&& indicatorFunction : _indicatorFunctions )
      f.add( indicatorFunction.a(), indicatorFunction.b(), -indicatorFunction.y() );

    return f;
  }

  /** Adds a scalar to all step function values */
  StepFunction operator+( I lambda ) const noexcept
  {
    StepFunction f;

    for( auto&& indicatorFunction : _indicatorFunctions )
      f.add( indicatorFunction.a(), indicatorFunction.b(), lambda + indicatorFunction.y() );

    return f;
  }

  /** Subtracts a scalar from all step function values */
  StepFunction operator-( I lambda ) const noexcept
  {
    return this->operator+( -lambda );
  }

  /** Multiplies the given step function with a scalar value */
  StepFunction operator*( I lambda ) const noexcept
  {
    StepFunction f;

    for( auto&& indicatorFunction : _indicatorFunctions )
      f.add( indicatorFunction.a(), indicatorFunction.b(), lambda * indicatorFunction.y() );

    return f;
  }

  /** Divides the given step function by a scalar value */
  StepFunction operator/( I lambda ) const
  {
    // TODO: What about division by zero?
    return this->operator*( 1/lambda );
  }

  /** Calculates the integral over the domain of the step function */
  I integral() const noexcept
  {
    if( _indicatorFunctions.empty() )
      return I();

    I value = I();

    for( auto&& f : _indicatorFunctions )
      value += f.integral();

    return value;
  }

  /** Calculates the unsigned integral raised to a certain power */
  I integral_p( I p ) const noexcept
  {
    if( _indicatorFunctions.empty() )
      return I();

    I value = I();

    for( auto&& f : _indicatorFunctions )
      value += f.integral_p( p );

    return std::pow( value, 1/p );
  }

  /** Calculates the absolute value of the function */
  StepFunction& abs() noexcept
  {
    std::set<IndicatorFunction> indicatorFunctions;

    for( auto&& f : _indicatorFunctions )
      indicatorFunctions.insert( IndicatorFunction( f.a(), f.b(), std::abs( f.y() ) ) );

    _indicatorFunctions.swap( indicatorFunctions );
    return *this;
  }

  template <class U, class V> friend std::ostream& operator<<( std::ostream&, const StepFunction<U, V>& f );

private:

  /** All indicator functions of the step function */
  std::set<IndicatorFunction> _indicatorFunctions;
};

// TODO: This does not need to be a friend function; it suffices to be
// implemented using the public interface of the class.
template <class D, class I> std::ostream& operator<<( std::ostream& o, const StepFunction<D, I>& f )
{
  for( auto&& indicatorFunction : f._indicatorFunctions )
  {
    o << indicatorFunction.a() << "\t" << indicatorFunction.y() << "\n"
      << indicatorFunction.b() << "\t" << indicatorFunction.y() << "\n";
  }

  return o;
}

} // namespace math

} // namespace aleph

#endif
