#ifndef ALEPH_MATH_STEP_FUNCTION_HH__
#define ALEPH_MATH_STEP_FUNCTION_HH__

#include <algorithm>
#include <iterator>
#include <ostream>
#include <set>
#include <stdexcept>

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

    bool operator<( const IndicatorFunction& other ) const
    {
      // Permit that intervals intersect in a single point, as this is
      // simplifies the composition of multiple step functions.
      return this->b() <= other.a();
    }

  private:
    T _a;
    T _b;
    T _y;
  };

  /**
    Auxiliary class for representing a point 'on' the step function;
    this assumes that the function consists of its non-zero points.
  */

  class Point
  {
  public:
    Point( T x, T y )
      : _x( x )
      , _y( y )
    {
    }

          T& x()       noexcept { return _x; }
    const T& x() const noexcept { return _x; }

          T& y()       noexcept { return _y; }
    const T& y() const noexcept { return _y; }

    bool operator<( const Point& other ) const
    {
      // Only compare by x coordinate; there must not be two or more
      // points with the same coordinate.
      return this->x() < other.x();
    }

  private:
    T _x = T();
    T _y = T();
  };

  /** Adds a new point to the step function */
  void add( T x, T y ) noexcept
  {
    _points.insert( Point(x,y) );
  }

  /** Adds a new indicator function to the step function */
  void add( T a, T b, T y ) noexcept
  {
    _indicatorFunctions.insert( IndicatorFunction(a,b,y) );
  }

  /** Returns the domain of the function */
  template <class OutputIterator> void domain( OutputIterator result )
  {
    for( auto&& f : _indicatorFunctions )
    {
      *result++ = f.a();
      *result++ = f.b();
    }
  }

  /** Returns the image of the function */
  template <class OutputIterator> void image( OutputIterator result )
  {
    for( auto&& f : _indicatorFunctions )
      *result++ = f.y();
  }

  /** Returns the function value at a certain position */
  T operator()( T x ) const noexcept
  {
    // TODO: Exploit the sorting order of the indicator functions and use the
    // nearest interval that contains the point.
    //
    // TODO: Handle *multiple* intervals...
    auto it = std::find_if( _indicatorFunctions.begin(), _indicatorFunctions.end(),
                            [&x] ( const IndicatorFunction& f )
                            {
                              return f.contains(x);
                            } );

    if( it != _indicatorFunctions.end() )
      return it->y();
    else
      return T();
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

    for( auto&& x : domain )
    {
      auto y1 = f(x);
      auto y2 = g(x);

      h.add( x, y1+y2 );
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
    if( _points.empty() )
      return T();

    auto cur = _points.begin();
    auto pre = _points.begin();

    std::advance( cur, 1 );

    T value = T();

    for( ; cur != _points.end(); )
    {
      auto c  = this->operator()( pre->x() );
      auto t  = cur->x() - pre->x();
      value  += c*t;

      pre = cur++;
    }

    return value;
  }

  template <class U> friend std::ostream& operator<<( std::ostream&, const StepFunction<U>& f );

private:

  /** All indicator functions of the step function */
  std::set<IndicatorFunction> _indicatorFunctions;

  /** All non-zero points of the step function */
  std::set<Point> _points;
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
