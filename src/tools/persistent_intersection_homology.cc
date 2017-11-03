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
#include <numeric>
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

template <class Functor> std::vector<DataType> extract( const PointCloud& pointCloud, Functor f )
{
  std::vector<DataType> values;
  values.reserve( pointCloud.size() );

  for( std::size_t i = 0; i < pointCloud.size(); i++ )
  {
    auto p = pointCloud[i];
    auto x = f( p.begin(), p.end() );

    values.push_back(x);
  }

  return values;
}

int main( int argc, char** argv )
{
  if( argc <= 2 )
    return -1;

  std::string inputPointCloud = argv[1];
  std::string inputCurvatures = "";
  DataType epsilon            = static_cast<DataType>( std::stod( argv[2] ) );

  if( argc >= 3 )
    inputCurvatures = argv[3];

  auto pointCloud = aleph::containers::load<DataType>( inputPointCloud );

  std::vector<DataType> singularityValues;
  singularityValues.reserve( pointCloud.size() );

  if( inputCurvatures.empty() == false )
  {
    auto curvatures   = aleph::containers::load<DataType>( inputCurvatures );
    singularityValues = extract( curvatures,
                                  [] ( auto begin, auto end )
                                  {
                                    return std::accumulate( begin, end, DataType() );
                                  } );
  }

  auto K
    = aleph::geometry::buildVietorisRipsComplex(
        NearestNeighbours( pointCloud ),
        epsilon,
        2 // FIXME: make configurable
  );

  std::cerr << "* Obtained Vietoris--Rips complex with " << K.size() << " simplices\n";

  decltype(K) K0, K1, K2, L;

  // Determine stratification ------------------------------------------
  //
  // There are two modes of operation here. First, if no singularity
  // values have been specified by the user, we employ the canonical
  // stratification based on skeletons.
  if( singularityValues.empty() )
  {
    std::cerr << "* Calculating skeletons...";

    K0 = aleph::topology::Skeleton()( 0, K );
    K1 = K0;
    K2 = K;

    std::cerr << "finished\n";

    std::cerr << "* Performing barycentric subdivision...";

    // Barycentric subdivision to ensure that the resulting complex is
    // flaglike in the sense of MacPherson et al.
    L = aleph::topology::BarycentricSubdivision()( K, [] ( std::size_t dimension ) { return dimension == 0 ? 0 : 0.5; } );

    {
      bool useMaximum                  = true;
      bool skipOneDimensionalSimplices = true;

      L.recalculateWeights( useMaximum, skipOneDimensionalSimplices );
      L.sort( aleph::topology::filtrations::Data<typename decltype(L)::ValueType>() ); // FIXME
    }

    std::cerr << "finished\n"
              << "* Subdivided simplicial complex has " << L.size() << " simplices\n";
  }

  // Else, we use the supplied singularity values to forbid parts of the
  // original data sets because they are too close to a singularity.
  else
  {
  }

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

  if( D1.size() >= 2 && D2.size() >= 2 )
  {
    std::ofstream out0( "/tmp/D_1_PH.txt" );
    std::ofstream out1( "/tmp/D_1_IH.txt" );

    D1[1].removeDiagonal();
    D2[1].removeDiagonal();

    out0 << D1[1] << "\n";
    out1 << D2[1] << "\n";
  }
}
