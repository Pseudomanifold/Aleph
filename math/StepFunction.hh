#ifndef ALEPH_MATH_STEP_FUNCTION_HH__
#define ALEPH_MATH_STEP_FUNCTION_HH__

#include <iterator>
#include <limits>
#include <ostream>
#include <set>
#include <stdexcept>

#include <cmath>

namespace aleph
{

namespace math
{

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

    I integral() const noexcept
    {
      return this->y() * static_cast<I>( ( this->b() - this->a() ) );
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

    auto prev = domain.begin();
    auto curr = domain.begin();

    I value = I();

    for( ; curr != domain.end(); )
    {
      if( prev != curr )
      {
        auto evaluationPoint = (*curr + *prev) / 2;

        auto y1 = f(evaluationPoint);
        auto y2 = g(evaluationPoint);

        if( y1+y2 != value )
        {
          h.add( *prev, *curr, y1+y2 );
          value = y1+y2;
        }
      }

      prev = curr++;
    }

    return h;
  }

  /** Multiplies the given step function with a scalar value */
  StepFunction operator*( I lambda ) const noexcept
  {
    auto f = *this;

    for( auto&& p : f._points )
      p.y() = p.y() * lambda;
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
