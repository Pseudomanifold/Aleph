
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

#include <iostream>
#include <string>

#include <getopt.h>

int main( int argc, char** argv )
{
  using DataType   = double;
  using PointCloud = aleph::containers::PointCloud<DataType>;
  using Distance   = aleph::distances::Euclidean<DataType>;

  #ifdef ALEPH_WITH_FLANN
    using FLANN = aleph::geometry::FLANN<PointCloud, Distance>;
  #endif

  static option commandLineOptions[] =
  {
    { "dimension"     , required_argument, nullptr, 'D' },
    { "descriptor"    , required_argument, nullptr, 'd' },
    { "epsilon"       , required_argument, nullptr, 'e' },
    { nullptr         , 0                , nullptr,  0  }
  };

  unsigned dimension     = 0;
  DataType epsilon       = DataType();
  std::string descriptor = std::string();

  {
    int option = 0;
    while( ( option = getopt_long( argc, argv, "D:d:e:", commandLineOptions, nullptr ) ) != -1 )
    {
      switch( option )
      {
      case 'D':
        dimension = static_cast<unsigned>( std::stoul( optarg ) );
        break;
      case 'd':
        descriptor = optarg;
        break;
      case 'e':
        epsilon = static_cast<DataType>( std::stod( optarg ) );
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

  auto dataDescriptorValues = aleph::estimateDensityDistanceToMeasure<Distance, PointCloud>( pointCloud, 10 );

  // TODO:
  //   - Data descriptor selection
  //   - Data descriptor modification according to some (?) criteria
  //   - Make unpaired simplex removal possible

  // Expansion ---------------------------------------------------------

  std::cerr << "* Expanding point cloud using epsilon=" << epsilon << "...";

  #ifdef ALEPH_WITH_FLANN
    FLANN flannWrapper( pointCloud );

    auto K
      = aleph::geometry::buildVietorisRipsComplex( flannWrapper,
                                                   epsilon,
                                                   unsigned( dimension ),
                                                   dataDescriptorValues.begin(), dataDescriptorValues.end() );

  #else
    aleph::geometry::BruteForce<PointCloud, Distance> flannWrapper( pointCloud );

    auto K
      = aleph::geometry::buildVietorisRipsComplex( flannWrapper,
                                                   epsilon,
                                                   unsigned( dimension ),
                                                   dataDescriptorValues.begin(), dataDescriptorValues.end() );
  #endif

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

    std::cout << "# Persistence diagram <" << input << ">\n"
              << "#\n"
              << "# Dimension: " << D.dimension() << "\n"
              << "# Entries  : " << D.size() << "\n"
              << D << "\n\n";
  }
}
