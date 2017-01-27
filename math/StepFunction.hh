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

template <class T> class StepFunction
{
public:

  /**
    Auxiliary class for representing an indicator function interval of the step
    function. Each indicator function is only non-zero within its interval, and
    zero outside.
  */

  class IndicatorFunction
  {
  public:
    IndicatorFunction( T a )
      : _a( a    )
      , _b( a    )
      , _y( T(1) )
    {
    }

    IndicatorFunction( T a, T b, T y )
      : _a( a )
      , _b( b )
      , _y( y )
    {
      if( a > b )
        throw std::runtime_error( "Invalid interval specified" );
    }

    const T& a() const noexcept { return _a; }
    const T& b() const noexcept { return _b; }

          T& y()       noexcept { return _y; }
    const T& y() const noexcept { return _y; }

    bool contains( T x ) const noexcept
    {
      return this->a() <= x && x <= this->b();
    }

    T integral() const noexcept
    {
      return this->y() * ( this->b() - this->a() );
    }

    T operator()( T x ) const noexcept
    {
      if( this->a() <= x && x <= this->b() )
        return this->y();
      else
        return T();
    }

    bool operator<( const IndicatorFunction& other ) const
    {
      // Permits that intervals intersect in a single point, as this is
      // simplifies the composition of multiple step functions.
      return this->b() <= other.a();
    }

  private:
    T _a;
    T _b;
    T _y;
  };

  /** Adds a new indicator function to the step function */
  void add( T a, T b, T y ) noexcept
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
  T operator()( T x ) const noexcept
  {
    T value = T();

    for( auto&& f : _indicatorFunctions )
    {
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

    std::set<T> domain;

    f.domain( std::inserter( domain, domain.begin() ) );
    g.domain( std::inserter( domain, domain.begin() ) );

    StepFunction<T> h;

    auto prev = domain.begin();
    auto curr = domain.begin();

    T value = T();

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
  StepFunction operator*( T lambda ) const noexcept
  {
    StepFunction<T> f = *this;

    for( auto&& p : f._points )
      p.y() = p.y() * lambda;
  }

  /** Divides the given step function by a scalar value */
  StepFunction operator/( T lambda ) const
  {
    // TODO: What about division by zero?
    return this->operator*( 1/lambda );
  }

  /** Calculates the integral over the domain of the step function */
  T integral() const noexcept
  {
    if( _indicatorFunctions.empty() )
      return T();

    T value = T();

    for( auto&& f : _indicatorFunctions )
      value += f.integral();

    return value;
  }

  template <class U> friend std::ostream& operator<<( std::ostream&, const StepFunction<U>& f );

private:

  /** All indicator functions of the step function */
  std::set<IndicatorFunction> _indicatorFunctions;
};

// TODO: This does not need to be a friend function; it suffices to be
// implemented using the public interface of the class.
template <class T> std::ostream& operator<<( std::ostream& o, const StepFunction<T>& f )
{
  for( auto&& p : f._points )
    o << p.x() << "\t" << p.y() << "\n";

  return o;
}

} // namespace math

} // namespace aleph

#endif
