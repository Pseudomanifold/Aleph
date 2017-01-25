#ifndef ALEPH_MATH_STEP_FUNCTION_HH__
#define ALEPH_MATH_STEP_FUNCTION_HH__

#include <algorithm>
#include <iterator>
#include <set>

namespace aleph
{

namespace math
{

template <class T> class StepFunction
{
public:

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

    T x() const noexcept { return _x; }
    T y() const noexcept { return _y; }

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

  /** Returns the function value at a certain position */
  T operator()( T x ) const noexcept
  {
    auto it = std::find_if( _points.begin(), _points.end(),
                            [&x] ( const Point& p )
                            {
                              return p.x() == x;
                            } );

    if( it != _points.end() )
      return it->y();
    else
      return T();
  }

  /** Calculates the integral over the domain of the step function */
  T integral() const noexcept
  {
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

private:

  /** All non-zero points of the step function */
  std::set<Point> _points;
};

} // namespace math

} // namespace aleph

#endif
