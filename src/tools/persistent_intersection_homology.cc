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

#include <aleph/math/Statistics.hh>

#include <aleph/persistenceDiagrams/PersistenceDiagram.hh>
#include <aleph/persistenceDiagrams/distances/Bottleneck.hh>

#include <aleph/persistentHomology/Calculation.hh>
#include <aleph/persistentHomology/PhiPersistence.hh>

#include <aleph/topology/BarycentricSubdivision.hh>
#include <aleph/topology/Filter.hh>
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

#include <getopt.h>

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

std::vector<DataType> standardizeValues( const std::vector<DataType>& data )
{
  auto mean   = aleph::math::sampleMean( data.begin(), data.end() );
  auto sdev   = aleph::math::sampleStandardDeviation( data.begin(), data.end() );
  auto result = data;

  std::transform( data.begin(), data.end(), result.begin(),
    [&mean, &sdev] ( DataType x )
    {
      return ( x - mean ) / sdev;
    }
  );

  return result;
}

int main( int argc, char** argv )
{
  DataType filterThreshold = DataType();
  bool invert              = false;
  bool standardize         = false;

  {
    static option commandLineOptions[] =
    {
      { "filter"     , required_argument, nullptr, 'f' },
      { "invert"     , no_argument      , nullptr, 'i' },
      { "standardize", no_argument      , nullptr, 's' },
      { nullptr      , 0                , nullptr,  0  }
    };

    int option = 0;
    while( ( option = getopt_long( argc, argv, "f:is", commandLineOptions, nullptr ) ) != -1 )
    {
      switch( option )
      {
      case 'f':
        filterThreshold = static_cast<DataType>( std::stod(optarg) );
        break;
      case 'i':
        invert = true;
        break;
      case 's':
        standardize = true;
        break;
      }
    }
  }

  if( ( argc - optind ) < 2 )
    return -1;

  std::string inputPointCloud = argv[optind++];
  std::string inputCurvatures = "";
  DataType epsilon            = static_cast<DataType>( std::stod( argv[optind++] ) );

  if( ( argc - optind ) >= 1 )
    inputCurvatures = argv[optind++];

  auto pointCloud = aleph::containers::load<DataType>( inputPointCloud );

  std::vector<DataType> singularityValues;
  singularityValues.reserve( pointCloud.size() );

  if( inputCurvatures.empty() == false )
  {
    std::cerr << "* Loading singularity values...";

    auto curvatures   = aleph::containers::load<DataType>( inputCurvatures );
    singularityValues = extract( curvatures,
                                  [] ( auto begin, auto end )
                                  {
                                    return std::accumulate( begin, end, DataType() );
                                  } );

    std::cerr << "finished\n";

    if( standardize )
    {
      std::cerr << "* Standardizing singularity values...";

      singularityValues = standardizeValues( singularityValues );

      std::cerr << "finished\n";
    }
  }

  unsigned dimension = 2;
  if( ( argc - optind ) >= 1 )
    dimension = static_cast<unsigned>( std::stoul( argv[optind++] ) );

  auto K
    = aleph::geometry::buildVietorisRipsComplex(
        NearestNeighbours( pointCloud ),
        epsilon,
        dimension
  );

  std::cerr << "* Obtained Vietoris--Rips complex with " << K.size() << " simplices\n";

  decltype(K) K0, K1, K2, K3, L;

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
    std::cerr << "* Using singularity values to filter complex...";

    aleph::topology::Filter filter;
    K0 = filter( K,
      [&filterThreshold,&invert,&singularityValues] ( auto s )
      {
        if( s.dimension() == 0 )
        {
          auto v = s[0];
          auto x = singularityValues.at(v);

          if( invert )
            return x > filterThreshold;
          else
            return x < filterThreshold;
        }

        return false;
      }
    );

    std::cerr << "finished\n"
              << "* Filtered 0-dimensional complex has " << K0.size() << " simplices\n";

    K1 = filter( K,
      [&filterThreshold,&invert,&singularityValues] ( auto s )
      {
        if( s.dimension() == 0 )
        {
          auto v = s[0];
          auto x = singularityValues.at(v);

          if( invert )
            return x > filterThreshold;
          else
            return x < filterThreshold;
        }
        else if( s.dimension() == 1 )
        {
          auto u = s[0];
          auto v = s[1];

          auto x = singularityValues.at(u);
          auto y = singularityValues.at(v);

          if( invert )
            return std::max(x,y) > filterThreshold;
          else
            return std::max(x,y) < filterThreshold;
        }

        return false;
      }
    );

    if( K.dimension() == 2 )
    {
      K1 = K0;
      K2 = K;
    }
    else
    {
      K2 = K1;
      K3 = K;
    }

    L  = K;
  }

  std::cerr << "* Calculating persistent homology...";

  auto D1 = aleph::calculatePersistenceDiagrams( K );

  std::cerr << "finished\n";

  std::cerr << "* Calculating intersection homology...";

  auto D2 = D1;
  if( K.dimension() == 2 )
    D2 = aleph::calculateIntersectionHomology( L, {K0,K1,K2}, aleph::PerversityGM( {0} ) );
  else
    D2 = aleph::calculateIntersectionHomology( L, {K0,K1,K2,K3}, aleph::PerversityGM( {0,1} ) );

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

  std::cerr << D1.size() << "," << D2.size() << "\n";

  if( D1.size() >= 3 && D2.size() >= 3 )
  {
    std::ofstream out0( "/tmp/D_2_PH.txt" );
    std::ofstream out1( "/tmp/D_2_IH.txt" );

    D1[2].removeDiagonal();
    D2[2].removeDiagonal();

    out0 << D1[2] << "\n";
    out1 << D2[2] << "\n";
  }
}
