#ifndef ALEPH_MATH_ALGEBRAIC_SPHERE_HH__
#define ALEPH_MATH_ALGEBRAIC_SPHERE_HH__

#include <algorithm>
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
    if( _coefficients.empty() || _coefficients.back() == T() )
      throw std::runtime_error( "Invalid coefficient set" );
  }

  /** Calculates and returns the centre of the sphere */
  std::vector<T> centre() const noexcept
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

    return c;
  }

  /** Calculates and returns the radius of the sphere */
  T radius() const noexcept
  {
    auto c = this->centre();
    auto n = std::inner_product( c.begin(), c.end(), c.begin(), T() );

    return std::sqrt( n - _coefficients.front() / _coefficients.back() );
  }

  /** Calculates and returns the Gaussian curvature of the sphere */
  T gaussianCurvature() const
  {
    return 1 / ( this->radius() * this->radius() );
  }

  /** Calculates and returns the mean curvature of the sphere */
  T meanCurvature() const
  {
    return 1 / this->radius();
  }

private:

  /** Sphere coefficients */
  std::vector<T> _coefficients;
};

} // namespace math

} // namespace aleph

#endif
