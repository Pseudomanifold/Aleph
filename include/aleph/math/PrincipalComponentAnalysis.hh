#ifndef ALEPH_MATH_PRINCIPAL_COMPONENT_ANALYSIS_HH__
#define ALEPH_MATH_PRINCIPAL_COMPONENT_ANALYSIS_HH__

#include <aleph/config/Eigen.hh>

#ifdef ALEPH_WITH_EIGEN
  #include <Eigen/Core>
  #include <Eigen/SVD>
#endif

#include <vector>

#include <cmath>

// These warnings can become a bit overzealous; the initialization done
// in the class is completely fine and will default to a struct that is
// properly initialized.
_Pragma( "GCC diagnostic push" )
_Pragma( "GCC diagnostic ignored \"-Wmissing-field-initializers\"" )

namespace aleph
{

namespace math
{

class PrincipalComponentAnalysis
{
public:

  template <class T > struct Result
  {
    std::vector< std::vector<T> > components;
    std::vector<T> singularValues;
  };

  // Main functor ------------------------------------------------------

  template <class T> Result<T> operator()( const std::vector< std::vector<T> >& data )
  {
#ifdef ALEPH_WITH_EIGEN
    if( data.empty() )
      return {};

    auto n = data.size();
    auto m = data.front().size();

    using Matrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
    using Vector = Eigen::Matrix<T, 1, Eigen::Dynamic>;
    using Index  = Eigen::Index;

    Matrix M(n,m);

    for( std::size_t row = 0; row < n; row++ )
      M.row( Index(row) ) = Vector::Map( &data[row][0], Index(m) );

    M  = M.rowwise() - M.colwise().mean();
    M /= std::sqrt( static_cast<T>( m ) );

    Eigen::JacobiSVD<Matrix> svd( M, Eigen::ComputeThinV );

    Result<T> result;

    {
      auto&& singularValues = svd.singularValues();
      result.singularValues.reserve( static_cast<std::size_t>( singularValues.size() ) );

      for( decltype( singularValues.size() ) i = 0; i < singularValues.size(); i++ )
        result.singularValues.push_back( singularValues( Index(i) ) );
    }

    {
      auto numSingularVectors = std::min( n, m );
      auto dimension          = m;

      result.components.resize( numSingularVectors,
                                std::vector<T>() );

      auto&& V = svd.matrixV();

      for( decltype(numSingularVectors) i = 0; i < numSingularVectors; i++ )
      {
        auto&& column = V.col( Index(i) );
        result.components[i].assign( column.data(), column.data() + dimension );
      }
    }

    return result;

#else
  // to quiet compiler warnings
  (void) data;
  return {};
#endif
  }
};

} // namespace math

} // namespace aleph

_Pragma( "GCC diagnostic pop" )

#endif
