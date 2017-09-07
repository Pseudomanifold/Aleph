#include <aleph/containers/PointCloud.hh>

#include <aleph/geometry/BetaSkeleton.hh>
#include <aleph/geometry/HeatKernel.hh>

#include <aleph/geometry/distances/Euclidean.hh>

#include <iostream>
#include <unordered_map>
#include <string>
#include <vector>

using DataType   = double;
using Distance   = aleph::distances::Euclidean<DataType>;
using PointCloud = aleph::containers::PointCloud<DataType>;

struct ScaleEstimationFunctor
{
  template <class SimplicialComplex> std::vector<DataType> operator()( const SimplicialComplex& K )
  {
    using Simplex    = typename SimplicialComplex::ValueType;
    using VertexType = typename Simplex::VertexType;
    using IndexType  = std::size_t;

    std::unordered_map<VertexType, IndexType> vertex_to_index;

    {
      std::vector<VertexType> vertices;
      K.vertices( std::back_inserter( vertices ) );

      IndexType index = IndexType();
      for( auto&& vertex : vertices )
        vertex_to_index[vertex] = index++;
    }

    std::vector<DataType> sumOfWeights( vertex_to_index.size() );
    std::vector<DataType> degree( vertex_to_index.size() );

    for( auto&& simplex : K )
    {
      if( simplex.dimension() != 1 )
        continue;

      auto&& u = simplex[0];
      auto&& v = simplex[1];
      auto&& w = simplex.data();
      auto&& i = vertex_to_index.at(u);
      auto&& j = vertex_to_index.at(v);

      degree[i] += 1;
      degree[j] += 1;

      sumOfWeights[i] += w;
      sumOfWeights[j] += w;
    }

    std::vector<DataType> scales( vertex_to_index.size() );

    for( std::size_t i = 0; i < vertex_to_index.size(); i++ )
      scales[i] = sumOfWeights[i] / degree[i];

    return scales;
  }
};

int main( int argc, char** argv )
{
  if( argc <= 1 )
    return -1;

  auto filename   = std::string( argv[1] );
  auto pointCloud = aleph::containers::load<DataType>( filename );

  std::cerr << "* Loaded point cloud with " << pointCloud.size() << " points\n";

  // Skeleton construction ---------------------------------------------

  // TODO: make configurable
  DataType beta = 1.0;

  std::cerr << "* Calculating beta-skeleton with beta = " << beta << "...";

  auto betaSkeleton
    = aleph::geometry::buildBetaSkeletonNaive( pointCloud,
                                               beta,
                                               Distance() );

  std::cerr << "...finished\n"
            << "* Simplical complex has " << betaSkeleton.size() << " simplices\n";

  // Scale estimation --------------------------------------------------

  ScaleEstimationFunctor sef;

  auto scalesBefore = sef( betaSkeleton );

  std::cerr << "* Initial scale information: ";
  for( auto&& s : scalesBefore )
    std::cerr << s << " ";
  std::cerr << "\n";
}
