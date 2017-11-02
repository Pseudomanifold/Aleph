#include <aleph/containers/DataDescriptors.hh>
#include <aleph/containers/DimensionalityEstimators.hh>
#include <aleph/containers/PointCloud.hh>

#include <aleph/config/FLANN.hh>

#ifdef ALEPH_WITH_FLANN
  #include <aleph/geometry/FLANN.hh>
#else
  #include <aleph/geometry/BruteForce.hh>
#endif

#include <aleph/geometry/SphereSampling.hh>
#include <aleph/geometry/VietorisRipsComplex.hh>

#include <aleph/geometry/distances/Euclidean.hh>

#include <aleph/persistenceDiagrams/PersistenceDiagram.hh>
#include <aleph/persistenceDiagrams/distances/Bottleneck.hh>

#include <aleph/persistentHomology/Calculation.hh>
#include <aleph/persistentHomology/PhiPersistence.hh>

#include <aleph/topology/BarycentricSubdivision.hh>
#include <aleph/topology/Simplex.hh>
#include <aleph/topology/SimplicialComplex.hh>
#include <aleph/topology/Skeleton.hh>

#include <aleph/topology/filtrations/Data.hh>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <set>
#include <vector>

using DataType           = double;
using VertexType         = unsigned;
using Distance           = aleph::geometry::distances::Euclidean<DataType>;
using PointCloud         = aleph::containers::PointCloud<DataType>;
using Simplex            = aleph::topology::Simplex<DataType, VertexType>;
using SimplicialComplex  = aleph::topology::SimplicialComplex<Simplex>;
using Filtration         = aleph::topology::filtrations::Data<Simplex>;
using PersistenceDiagram = aleph::PersistenceDiagram<DataType>;

#ifdef ALEPH_WITH_FLANN
  using NearestNeighbours = aleph::geometry::FLANN<PointCloud, Distance>;
#else
  using NearestNeighbours = aleph::geometry::BruteForce<PointCloud, Distance>;
#endif

int main( int argc, char** argv )
{
  if( argc <= 2 )
    return -1;

  std::string inputPointCloud = argv[1];
  std::string inputCurvatures = argv[2];

  auto pointCloud = aleph::containers::load<DataType>( inputPointCloud );
  auto curvatures = aleph::containers::load<DataType>( inputCurvatures );

  auto K
    = aleph::geometry::buildVietorisRipsComplex(
        NearestNeighbours( pointCloud ),
        DataType( 0.50 ), // FIXME: make configurable
        2                 // FIXME: make configurable
  );

  std::cerr << "* Obtained Vietoris--Rips complex with " << K.size() << " simplices\n";

  std::cerr << "* Calculating skeletons...";

  auto K0 = aleph::topology::Skeleton()( 0, K );
  auto K1 = K0;
  auto K2 = K;

  std::cerr << "finished\n";

  std::cerr << "* Performing barycentric subdivision...";

  // Barycentric subdivision to ensure that the resulting complex is
  // flaglike in the sense of MacPherson et al.
  auto L
    = aleph::topology::BarycentricSubdivision()( K, [] ( std::size_t dimension ) { return dimension == 0 ? 0 : 0.5; } );

  {
    bool useMaximum                  = true;
    bool skipOneDimensionalSimplices = true;

    L.recalculateWeights( useMaximum, skipOneDimensionalSimplices );
    L.sort( aleph::topology::filtrations::Data<typename decltype(L)::ValueType>() ); // FIXME
  }

  std::cerr << "finished\n"
            << "* Subdivided simplicial complex has " << L.size() << " simplices\n";

  bool useOriginalIndexing = true;

  std::cerr << "* Calculating persistent homology...";

  auto D1 = aleph::calculatePersistenceDiagrams( K );

  std::cerr << "finished\n";

  std::cerr << "* Calculating intersection homology...";

  auto D2 = aleph::calculateIntersectionHomology( L, {K0,K1,K2}, aleph::PerversityGM( {0} ), useOriginalIndexing );

  std::cerr << "finished\n";

  {
    std::ofstream out0( "/tmp/D_0_PH.txt" );
    std::ofstream out1( "/tmp/D_0_IH.txt" );

    D1.front().removeDiagonal();
    D2.front().removeDiagonal();

    out0 << D1.front() << "\n";
    out1 << D2.front() << "\n";
  }

  {
    std::ofstream out0( "/tmp/D_1_PH.txt" );
    std::ofstream out1( "/tmp/D_1_IH.txt" );

    D1.back().removeDiagonal();
    D2.back().removeDiagonal();

    out0 << D1.back() << "\n";
    out1 << D2.back() << "\n";
  }
}
