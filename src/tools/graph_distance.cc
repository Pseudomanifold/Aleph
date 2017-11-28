#include <aleph/geometry/HeatKernel.hh>

#include <aleph/topology/io/GML.hh>

#include <aleph/topology/Simplex.hh>
#include <aleph/topology/SimplicialComplex.hh>

#include <aleph/utilities/Filesystem.hh>

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

    std::sort( spectrum1.begin(), spectrum1.end() );
    std::sort( spectrum2.begin(), spectrum2.end() );

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

  std::vector<std::string> filenames;

  for( int i = 1; i < argc; i++ )
    filenames.push_back( argv[i] );

  aleph::topology::io::GMLReader reader;

  for( auto&& filename : filenames )
  {
    std::cerr << "* Processing '" << filename << "'...";

    SimplicialComplex K;
    reader( filename, K );

    K.sort();

    // TODO:
    //  - Extract vertex data
    //  - Use vertex data in Laplacian matrix

    auto L = aleph::geometry::weightedLaplacianMatrix( K );

    std::cerr << "finished\n";
  }
}
