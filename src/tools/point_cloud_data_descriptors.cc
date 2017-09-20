
/*
  This is a tool shipped by 'Aleph - A Library for Exploring Persistent
  Homology'.

  It closely follows the paper

    Persistent Homology for the Evaluation of Dimensionality Reduction Schemes
    Bastian Rieck, Heike Leitte
    Computer Graphics Forum, Volume 34, Issue 3, pp. 431--440

  and implements a multitude of data descriptors that may be used during
  the expansion of a point cloud.

  The application knows two modes.

  1. Calculation of data descriptors and persistent homology
  2. Calculation of persistent homology based on existing data

  As of now, only the first mode is implemented.
*/

#include <aleph/config/FLANN.hh>

#include <aleph/containers/DataDescriptors.hh>
#include <aleph/containers/PointCloud.hh>

#include <aleph/geometry/BruteForce.hh>
#include <aleph/geometry/FLANN.hh>
#include <aleph/geometry/VietorisRipsComplex.hh>

#include <aleph/geometry/distances/Euclidean.hh>

#include <aleph/persistenceDiagrams/Calculation.hh>
#include <aleph/persistenceDiagrams/PersistenceDiagram.hh>

#include <aleph/persistenceDiagrams/io/JSON.hh>

#include <aleph/persistentHomology/Calculation.hh>

#include <algorithm>
#include <iostream>
#include <string>

#include <getopt.h>

using DataType   = double;
using PointCloud = aleph::containers::PointCloud<DataType>;
using Distance   = aleph::geometry::distances::Euclidean<DataType>;

#ifdef ALEPH_WITH_FLANN
  using Wrapper = aleph::geometry::FLANN<PointCloud, Distance>;
#else
  using Wrapper = aleph::geometry::BruteForce<PointCloud, Distance>;
#endif

void normalizeValues( std::vector<DataType>& values )
{
  if( values.empty() )
    return;

  auto minmax = std::minmax_element( values.begin(), values.end() );
  auto min    = *minmax.first;
  auto max    = *minmax.second;

  if( min == max )
    return;

  auto range = max - min;
  for( auto&& value : values )
    value = (value - min) / range;
}

void invertValues( std::vector<DataType>& values )
{
  if( values.empty() )
    return;

  auto max = *std::max_element( values.begin(), values.end() );
  for( auto&& value : values )
    value = max - value;
}

std::vector<DataType> calculateDataDescriptor( const std::string& name, const PointCloud& pointCloud, unsigned k, double h, unsigned p )
{
  if( name == "density" )
  {
    #ifdef ALEPH_WITH_FLANN
      return aleph::containers::estimateDensityDistanceToMeasure<Distance, PointCloud, Wrapper>( pointCloud, k );
    #else
      // The function will automatically fall back to the default
      // wrapper.
      return aleph::containers::estimateDensityDistanceToMeasure<Distance>( pointCloud, k );
    #endif
  }
  else if( name == "eccentricity" )
    return aleph::containers::eccentricities<Distance>( pointCloud, p );
  else if( name == "gaussian" )
    return aleph::containers::estimateDensityTruncatedGaussian( pointCloud, h );

  return {};
}

void usage()
{
  std::cerr << "Usage: point_cloud_data_descriptors [--bandwidth=H] [--dimension=D]\n"
            << "                                    [--descriptor=DESC]\n"
            << "                                    [--epsilon=EPS] [--k=k]\n"
            << "                                    [--invert] [--normalize]\n"
            << "                                    [--power=p]\n"
            << "                                    [--remove-unpaired] FILENAME\n"
            << "\n"
            << "Performs Vietoris--Rips expansion on the specified point cloud and\n"
            << "calculates its persistent homology based on the values of one data\n"
            << "descriptor. The expansion process uses an epsilon value of EPS and\n"
            << "a maximum dimension of D\n"
            << "\n"
            << "The following data descriptors are available as a name for DESC:\n"
            << "- density: uses distance to a measure density estimation. Notice\n"
            << "           that this descriptor queries the k nearest neighbours\n"
            << "           of a data point. By default, k=10, but this behaviour\n"
            << "           can be changed.\n"
            << "\n"
            << "- eccentricity: calculates eccentricity values for every point;\n"
            << "                the eccentricity measures the centrality of all\n"
            << "                points in the point cloud. Every value is taken\n"
            << "                to the p-th power, with p=2 by default. Specify\n"
            << "                p=0 in order to calculate maximum eccentricity.\n"
            << "\n"
            << "- gaussian: uses a truncated Gaussian density estimator with a\n"
            << "            bandwidth of h. By default, h=0.01.\n"
            << "\n"
            << "Several flags permit some control over the calculations:\n"
            << "--invert: inverts data descriptor values. This is useful for the\n"
            << "          eccentricity descriptor, for example, because it uses\n"
            << "          small values to indicate very central points.\n"
            << "\n"
            << "--normalize: normalizes data descriptor values to [0,1]\n"
            << "\n"
            << "--remove-unpaired: removes all unpaired simplices, thereby making\n"
            << "                   sure that all features have finite persistence\n"
            << "                   values\n"
            << "\n"
            << "Abbreviations of the command-line arguments specified above\n"
            << "are also supported:\n"
            << "  -b: bandwidth\n"
            << "  -D: dimension\n"
            << "  -d: descriptor\n"
            << "  -e: epsilon\n"
            << "  -k: number of nearest neighbours\n"
            << "  -i: invert values (no argument)\n"
            << "  -n: normalize values (no argument)\n"
            << "  -p: power for eccentricity calculation\n"
            << "  -r: remove unpaired simplices (no argument)\n"
            << "\n";
}

int main( int argc, char** argv )
{
  static option commandLineOptions[] =
  {
    { "bandwidth"      , required_argument, nullptr, 'b' },
    { "dimension"      , required_argument, nullptr, 'D' },
    { "descriptor"     , required_argument, nullptr, 'd' },
    { "epsilon"        , required_argument, nullptr, 'e' },
    { "k"              , required_argument, nullptr, 'k' },
    { "invert"         , no_argument      , nullptr, 'i' },
    { "normalize"      , no_argument      , nullptr, 'n' },
    { "power"          , required_argument, nullptr, 'p' },
    { "remove-unpaired", no_argument      , nullptr, 'r' },
    { nullptr          , 0                , nullptr,  0  }
  };

  unsigned dimension     = 0;           // default dimension (point cloud expansion)
  double h               = 0.01;        // default bandwidth (Gaussian estimator)
  unsigned k             = 10;          // default number of neighbours (density estimator)
  unsigned p             = 2;           // default power (eccentricity estimator)
  DataType epsilon       = DataType();  // default epsilon (point cloud expansion)
  std::string descriptor = "density";   // default data descriptor

  bool normalizeDataDescriptorValues = false;
  bool invertDataDescriptorValues    = false;
  bool removeUnpairedSimplices       = false;

  {
    int option = 0;
    while( ( option = getopt_long( argc, argv, "b:D:d:e:k:inp:r", commandLineOptions, nullptr ) ) != -1 )
    {
      switch( option )
      {
      case 'b':
        h = std::stod( optarg );
        break;
      case 'D':
        dimension = static_cast<unsigned>( std::stoul( optarg ) );
        break;
      case 'd':
        descriptor = optarg;
        break;
      case 'e':
        epsilon = static_cast<DataType>( std::stod( optarg ) );
        break;
      case 'k':
        k = static_cast<unsigned>( std::stoul( optarg ) );
        break;
      case 'i':
        invertDataDescriptorValues = true;
        break;
      case 'n':
        normalizeDataDescriptorValues = true;
        break;
      case 'p':
        p = static_cast<unsigned>( std::stoul( optarg ) );
        break;
      case 'r':
        removeUnpairedSimplices = true;
        break;
      }
    }
  }

  if( argc - optind <= 0 )
  {
    usage();
    return -1;
  }

  std::string input = argv[optind++];
  auto pointCloud   = aleph::containers::load<DataType>( input );

  if( dimension == 0 )
    dimension = static_cast<unsigned>( pointCloud.dimension() ) + 1;

  std::cerr << "* Obtained point cloud of dimension " << pointCloud.dimension() << " with " << pointCloud.size() << " points\n";

  auto dataDescriptorValues
    = calculateDataDescriptor( descriptor,
                               pointCloud,
                               k,
                               h,
                               p );

  if( invertDataDescriptorValues )
    invertValues( dataDescriptorValues );

  if( normalizeDataDescriptorValues )
    normalizeValues( dataDescriptorValues );

  // Expansion ---------------------------------------------------------

  std::cerr << "* Expanding point cloud using epsilon=" << epsilon << "...";

  Wrapper wrapper( pointCloud );

  auto K
    = aleph::geometry::buildVietorisRipsComplex( wrapper,
                                                 epsilon,
                                                 unsigned( dimension ),
                                                 dataDescriptorValues.begin(), dataDescriptorValues.end() );

  std::cerr << "finished\n"
            << "* Expanded simplicial complex has " << K.size() << " simplices\n";

  // Persistence diagram calculation ------------------------------------

  auto diagrams
    = aleph::calculatePersistenceDiagrams( K );

  std::cerr << "finished\n"
            << "* Obtained " << diagrams.size() << " persistence diagrams\n";

  std::cout << "{\n"
            << "\"diagrams\": " << "[\n";

  for( auto it = diagrams.begin(); it != diagrams.end(); ++it )
  {
    if( it != diagrams.begin() )
      std::cout << ",\n";

    auto&& D = *it;

    D.removeDiagonal();

    if( removeUnpairedSimplices )
      D.removeUnpaired();

    aleph::io::writeJSON( std::cout, D, input );
  }

  std::cout << "\n"
            << "  ]\n"
            << "}\n";
}
