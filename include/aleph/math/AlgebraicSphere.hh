#ifndef ALEPH_MATH_ALGEBRAIC_SPHERE_HH__
#define ALEPH_MATH_ALGEBRAIC_SPHERE_HH__

#include <algorithm>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <vector>

#include <cmath>

namespace aleph
{

namespace math
{

/**
  @class AlgebraicSphere
  @brief Models an algebraic sphere from a set of coefficients

  This class follows the concepts defined in the paper *Direct
  Least-Squares Fitting of Algebraic Surfaces* by V. Pratt and
  permits the calculation of some relevant properties, such as
  the *curvature* of the sphere.

  @see http://boole.stanford.edu/pub/fit.pdf
*/

template <class T> class AlgebraicSphere
{
public:

  /** Creates a new algebraic sphere from a set of coefficients */
  template <class InputIterator> AlgebraicSphere( InputIterator begin, InputIterator end )
    : _coefficients( begin, end )
  {
    if( _coefficients.size() < 3 )
      throw std::runtime_error( "At least three coefficients are required" );
  }

  /** Calculates and returns the centre of the sphere */
  std::vector<T> centre() const noexcept
  {
    if( _centre.empty() )
    {
      auto s = _coefficients.back();
      auto c = std::vector<T>( _coefficients.begin() + 1,
                               _coefficients.end()   - 1);

      std::transform( c.begin(), c.end(), c.begin(),
        [&s] ( T x )
        {
          return x / (2*s);
        }
      );

      _centre = c;
    }

    return _centre;
  }

  /** Calculates and returns the radius of the sphere */
  T radius() const noexcept
  {
    if( !std::isfinite( _radius ) )
    {
      auto c  = this->centre();
      auto n  = std::inner_product( c.begin(), c.end(), c.begin(), T() );
      _radius = std::sqrt( n - _coefficients.front() / _coefficients.back() );
    }

    return _radius;
  }

  /** Calculates and returns the Gaussian curvature of the sphere */
  T gaussianCurvature() const
  {
    auto r = this->radius();

    // Gracefully handle degenerate cases for which the sphere
    // degenerates into a plane.
    if( r > 0 )
      return 1 / ( r*r );
    else
      return T();
  }

  /** Calculates and returns the mean curvature of the sphere */
  T meanCurvature() const
  {
    auto r = this->radius();

    // Gracefully handle degenerate cases for which the sphere
    // degenerates into a plane.
    if( r > 0 )
      return 1 / r;
    else
      return T();
  }

private:

  /** Radius */
  mutable T _radius = std::numeric_limits<T>::quiet_NaN();

  /** Centre */
  mutable std::vector<T> _centre;

  /** Sphere coefficients */
  std::vector<T> _coefficients;
};

} // namespace math

} // namespace aleph

#endif
