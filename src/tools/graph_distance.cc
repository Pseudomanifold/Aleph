#include <aleph/config/Eigen.hh>

#include <aleph/geometry/HeatKernel.hh>

#include <aleph/topology/io/GML.hh>

#include <aleph/topology/Simplex.hh>
#include <aleph/topology/SimplicialComplex.hh>

#include <aleph/utilities/Filesystem.hh>

#ifdef ALEPH_WITH_EIGEN
  #include <Eigen/Eigenvalues>
#endif

#include <algorithm>
#include <iostream>
#include <map>
#include <numeric>
#include <regex>
#include <string>
#include <vector>

// Auxiliary class for storing the spectrum of a graph, i.e. the set of
// eigenvalues. This class offers a simple distance calculation.
template <class T> class Spectrum
{
public:

  template <class InputIterator> Spectrum( InputIterator begin, InputIterator end )
    : _eigenvalues( begin, end )
  {
  }

  std::size_t size() const noexcept
  {
    return _eigenvalues.size();
  }

  T distance( const Spectrum& other ) const noexcept
  {
    auto n         = this->size();
    auto m         = other.size();
    auto spectrum1 = this->_eigenvalues;
    auto spectrum2 = other._eigenvalues;

    spectrum1.resize( std::max(n,m) );
    spectrum2.resize( std::max(n,m) );

    std::sort( spectrum1.begin(), spectrum1.end(), std::greater<T>() );
    std::sort( spectrum2.begin(), spectrum2.end(), std::greater<T>() );

    // Squared Euclidean distance
    return std::inner_product( spectrum1.begin(), spectrum1.end(),
                               spectrum2.begin(),
                               T(0),
                               std::plus<T>(),
                               []( T x,T y )
                               {
                                 return ( x-y )*( x-y );
                               }
    );
  }

private:
  std::vector<T> _eigenvalues;
};

int main( int argc, char** argv )
{
  if( argc <= 1 )
    return -1;

  using DataType          = double;
  using VertexType        = unsigned short;
  using Simplex           = aleph::topology::Simplex<DataType, VertexType>;
  using SimplicialComplex = aleph::topology::SimplicialComplex<Simplex>;
  using Spectrum          = Spectrum<DataType>;

  std::vector<std::string> filenames;

  for( int i = 1; i < argc; i++ )
    filenames.push_back( argv[i] );

  aleph::topology::io::GMLReader reader;

  std::vector<Spectrum> spectra;
  spectra.reserve( filenames.size() );

  for( auto&& filename : filenames )
  {
    std::cerr << "* Processing '" << filename << "'...";

    SimplicialComplex K;
    reader( filename, K );

    K.sort();

    // TODO:
    //  - Extract vertex data
    //  - Use vertex data in Laplacian matrix

    auto L           = aleph::geometry::weightedLaplacianMatrix( K );
    using MatrixType = decltype(L);

    Eigen::SelfAdjointEigenSolver<MatrixType> solver;
    solver.compute( L );

    auto eigenvalues = solver.eigenvalues();

    {
      std::vector<DataType> data;
      for( decltype( eigenvalues.size() ) i = 0; i < eigenvalues.size(); i++ )
        data.push_back( eigenvalues( i ) );

      spectra.push_back( Spectrum( data.begin(), data.end() ) );
    }

    std::cerr << "finished\n";
  }

  // Calculate spectral distances --------------------------------------

  std::vector< std::vector<DataType> > distances( spectra.size(),
    std::vector<DataType>( spectra.size() ) );

  for( std::size_t i = 0; i < spectra.size(); i++ )
  {
    for( std::size_t j = i+1; j < spectra.size(); j++ )
    {
      distances[i][j] = spectra[i].distance( spectra[j] );
      distances[j][i] = distances[i][j];
    }
  }

  for( std::size_t i = 0; i < spectra.size(); i++ )
  {
    for( std::size_t j = 0; j < spectra.size(); j++ )
    {
      if( j != 0 )
        std::cout << " ";

      std::cout << distances[i][j];
    }

    std::cout << "\n";
  }

  std::cout << "\n\n";
}
