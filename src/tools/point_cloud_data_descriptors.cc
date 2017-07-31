
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

#include <aleph/persistentHomology/Calculation.hh>

#include <algorithm>
#include <iostream>
#include <string>

#include <getopt.h>

using DataType   = double;
using PointCloud = aleph::containers::PointCloud<DataType>;
using Distance   = aleph::distances::Euclidean<DataType>;

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

std::vector<DataType> calculateDataDescriptor( const std::string& name, const PointCloud& pointCloud, unsigned k, double h )
{
  if( name == "density" )
  {
    #ifdef ALEPH_WITH_FLANN
      return aleph::estimateDensityDistanceToMeasure<Distance, PointCloud, Wrapper>( pointCloud, k );
    #else
      // The function will automatically fall back to the default
      // wrapper.
      return aleph::estimateDensityDistanceToMeasure<Distance>( pointCloud, k );
    #endif
  }
  else if( name == "eccentricity" )
    return aleph::eccentricities<Distance>( pointCloud, k );
  else if( name == "gaussian" )
    return aleph::estimateDensityTruncatedGaussian( pointCloud, h );

  return {};
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
    { "remove-unpaired", no_argument      , nullptr, 'r' },
    { nullptr          , 0                , nullptr,  0  }
  };

  unsigned dimension     = 0;           // default dimension (point cloud expansion)
  double h               = 0.01;        // default bandwidth (Gaussian estimator)
  unsigned k             = 10;          // default number of neighbours (density estimator)
  DataType epsilon       = DataType();  // default epsilon (point cloud expansion)
  std::string descriptor = "density";   // default data descriptor

  bool normalizeDataDescriptorValues = false;
  bool invertDataDescriptorValues    = false;
  bool removeUnpairedSimplices       = false;

  {
    int option = 0;
    while( ( option = getopt_long( argc, argv, "b:D:d:e:k:inr", commandLineOptions, nullptr ) ) != -1 )
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
      case 'r':
        removeUnpairedSimplices = true;
        break;
      }
    }
  }

  if( argc - optind <= 0 )
    return -1;

  std::string input = argv[optind++];
  auto pointCloud   = aleph::containers::load<DataType>( input );

  if( dimension == 0 )
    dimension = static_cast<unsigned>( pointCloud.dimension() ) + 1;

  std::cerr << "* Obtained point cloud of dimension " << pointCloud.dimension() << " with " << pointCloud.size() << " points\n";

  auto dataDescriptorValues
    = calculateDataDescriptor( descriptor,
                               pointCloud,
                               k,
                               h );

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

  for( auto&& D : diagrams )
  {
    D.removeDiagonal();

    if( removeUnpairedSimplices )
      D.removeUnpaired();

    std::cout << "# Persistence diagram <" << input << ">\n"
              << "#\n"
              << "# Dimension: " << D.dimension() << "\n"
              << "# Entries  : " << D.size() << "\n"
              << D << "\n\n";
  }
}
