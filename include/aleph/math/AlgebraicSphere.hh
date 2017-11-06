#ifndef ALEPH_MATH_ALGEBRAIC_SPHERE_HH__
#define ALEPH_MATH_ALGEBRAIC_SPHERE_HH__

#include <algorithm>
#include <numeric>
#include <vector>

#include <cmath>

namespace aleph
{

namespace math
{

template <class T> class AlgebraicSphere
{
public:
  template <class InputIterator> AlgebraicSphere( InputIterator begin, InputIterator end )
    : _coefficients( begin, end )
  {
  }

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

  T radius() const noexcept
  {
    auto c = this->centre();
    auto n = std::inner_product( c.begin(), c.end(), c.begin(), T() );

    return std::sqrt( n - _coefficients.front() / _coefficients.back() );
  }

  T gaussianCurvature() const
  {
    return 1 / ( this->radius() * this->radius() );
  }

  T meanCurvature() const
  {
    return 1 / this->radius();
  }

private:
  std::vector<T> _coefficients;
};

} // namespace math

} // namespace aleph

#endif
