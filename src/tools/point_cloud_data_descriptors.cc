
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

#include <aleph/distances/Euclidean.hh>

#include <aleph/geometry/BruteForce.hh>
#include <aleph/geometry/FLANN.hh>
#include <aleph/geometry/VietorisRipsComplex.hh>

#include <aleph/persistenceDiagrams/Calculation.hh>
#include <aleph/persistenceDiagrams/PersistenceDiagram.hh>

#include <aleph/persistentHomology/Calculation.hh>

#include <iostream>
#include <string>

int main( int argc, char** argv )
{
  if( argc <= 2 )
    return -1;

  using DataType   = double;
  using PointCloud = aleph::containers::PointCloud<DataType>;
  using Distance   = aleph::distances::Euclidean<DataType>;
  using FLANN      = aleph::geometry::FLANN<PointCloud, Distance>;

  std::string input = argv[1];
  auto pointCloud   = aleph::containers::load<DataType>( input );
  auto dimension    = pointCloud.dimension() + 1;
  auto epsilon      = DataType();

  auto dataDescriptorValues = aleph::estimateDensityDistanceToMeasure<Distance, PointCloud>( pointCloud, 10 );

  // TODO:
  //   - Data descriptor selection
  //   - Data descriptor modification according to some (?) criteria
  //   - Make unpaired simplex removal possible

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
