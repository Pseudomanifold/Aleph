#ifndef ALEPH_MATH_STEP_FUNCTION_HH__
#define ALEPH_MATH_STEP_FUNCTION_HH__

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

private:

  /** All non-zero points of the step function */
  std::set<Point> _points;
};

} // namespace math

} // namespace aleph

#endif
