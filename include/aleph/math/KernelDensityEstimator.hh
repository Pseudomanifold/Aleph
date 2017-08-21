#ifndef ALEPH_MATH_KERNEL_DENSITY_ESTIMATOR_HH__
#define ALEPH_MATH_KERNEL_DENSITY_ESTIMATOR_HH__

#include <functional>
#include <iterator>
#include <numeric>
#include <type_traits>

#include <cmath>

namespace aleph
{

namespace math
{

/**
  @namespace kernels
  @brief     Kernels for the kernel density estimator class
*/

namespace kernels
{

/**
  @class Gaussian
  @brief Simple Gaussian kernel
*/

class Gaussian
{
public:
  Gaussian() = default;

  Gaussian( double sigma )
    : _sigma( sigma )
  {
  }

  template <class T> double operator()( T value ) const
  {
    return 1.0 / ( std::sqrt( 2.0 * M_PI ) * _sigma ) * std::exp( -0.5 * ( value * value ) / ( _sigma * _sigma ) );
  }

private:
  double _sigma = 1.0;
};

/**
  @class Epanechnikov
  @brief Simple Epanechnikov kernel
*/

class Epanechnikov
{
public:
  template <class T> double operator()( T value ) const
  {
    return std::abs( value ) <= static_cast<T>( 1.0 ) ? 0.75 * ( 1.0 - value*value ) : 0.0;
  }
};

} // namespace kernels


/**
  @namespace norms
  @brief     Norms for kernel density estimation
*/

namespace norms
{

/**
  @class Identity
  @brief Identity norm; returns a value unmodified
*/

class Identity
{
public:
  template <class T> T operator()( T value ) const
  {
    return value;
  }
};

/**
  @class Euclidean
  @brief Euclidean norm
*/

class Euclidean
{
public:
  template <class T> T operator()( T value ) const
  {
    auto squaredNorm
      = std::inner_product( std::begin( value ), std::end( value ),
                            std::begin( value ),
                            T() );

    return std::sqrt( squaredNorm );
  }
};

} // namespace norms

/**
  @class KernelDensityEstimator
  @brief Kernel density estimator class

  This class describes a generic kernel density estimator that works for
  univariate and multivariate data. The estimator is highly-configurable
  and permits the following settings:

  - Kernel selection

  - Norm selection (only relevant for multivariate KDE)

  - Difference selection in order to specify how differences between data
    points are being calculated
*/

class KernelDensityEstimator
{
public:

  /**
    Creates a new kernel density estimator with a given bandwidth that
    is capable of handling data of a certain dimensionality. Note that
    the dimensionality parameter is only used to scale the results. It
    has no bearing on the actual calculation.

    @param bandwidth Bandwidth
    @param dimension Dimension
  */

  KernelDensityEstimator( double bandwidth, unsigned dimension )
    : _bandwidth( bandwidth )
    , _dimension( dimension )
  {
  }

  /**
    Evaluates the kernel density estimator at a given position, using
    a certain set of parameters. The evaluator function is written so
    as to make type deduction easy for the compiler. For the 1D case,
    everything works automatically, provided a suitable kernel exists
    for the given data type.

    @param begin      Input iterator to begin of data range
    @param end        Input iterator to end of data range

    @param x          Point at which to evaluate the kernel density estimator

    @param kernel     Kernel to use for calculating density estimates

    @param norm       Norm to use for multivariate estimates; in the univariate
                      case, this parameter has no effect.

    @param difference Functor that specifies how differences between input data
                      are calculated.

    @tparam InputIterator Input iterator
    @tparam DataType      Type of the input data, i.e. `double` or `std::vector<double>`
    @tparam Kernel        Kernel to use for calculating density estimates
    @tparam Norm          Norm to use for multivariate estimates
    @tparam Difference    Functor for calculating differences for the input data type
  */

  template <
    class InputIterator,
    class DataType,
    class Kernel,
    class Norm       = norms::Identity,
    class Difference = std::minus<DataType>
  > double operator()( InputIterator begin, InputIterator end,
                       DataType x,
                       Kernel kernel = Kernel(),
                       Norm norm = Norm(),
                       Difference difference = Difference() ) const
  {
    static_assert( std::is_same<typename std::iterator_traits<InputIterator>::value_type, DataType>::value,
                   "Input iterator value type and data type must match" );

    double value = 0.0;

    for( InputIterator it = begin; it != end; ++it )
    {
      double norm_   = norm( difference( *it, x ) );
      double kernel_ = kernel( norm_ / _bandwidth );

      value += kernel_;
    }

    value /= std::pow( _bandwidth, static_cast<double>( _dimension ) );
    value /= static_cast<double>( std::distance( begin, end ) );

    return value;
  }

private:
  double   _bandwidth;
  unsigned _dimension;
};

} // namespace math

} // namespace aleph

#endif
