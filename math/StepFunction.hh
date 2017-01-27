#ifndef ALEPH_MATH_STEP_FUNCTION_HH__
#define ALEPH_MATH_STEP_FUNCTION_HH__

#include <algorithm>
#include <iterator>
#include <ostream>
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

  /** Returns the domain of the function */
  template <class OutputIterator> void domain( OutputIterator result )
  {
    for( auto&& p : _points )
      *result++ = p.x();
  }

  /** Returns the image of the function */
  template <class OutputIterator> void image( OutputIterator result )
  {
    for( auto&& p : _points )
      *result++ = p.y();
  }

  /** Returns the function value at a certain position */
  T operator()( T x ) const noexcept
  {
    // Find the point that is nearest to the query point and use its value as
    // the result. I don't want to do any interpolation here.
    auto it = std::lower_bound( _points.begin(), _points.end(),
                                Point( x, T() ) );

    if( it != _points.end() )
      return it->y();
    else
      return T();
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
