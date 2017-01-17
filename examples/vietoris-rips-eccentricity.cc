#include "containers/DataDescriptors.hh"
#include "containers/PointCloud.hh"

#include "distances/Euclidean.hh"

#include "geometry/FLANN.hh"
#include "geometry/VietorisRipsComplex.hh"

#include "utilities/String.hh"

#include "persistentHomology/Calculation.hh"

#include <iostream>
#include <string>

using namespace aleph;
using namespace distances;
using namespace geometry;
using namespace utilities;

int main( int argc, char** argv )
{
  if( argc <= 2 )
    return -1;

  using DataType   = double;
  using PointCloud = PointCloud<DataType>;
  using Distance   = Euclidean<DataType>;

  std::string input = argv[1];
  auto pointCloud   = load<DataType>( input );
  auto dimension    = pointCloud.dimension() + 1;
  auto epsilon      = convert<DataType>( argv[2] );
  auto order        = 1u;

  if( argc >= 4 )
    dimension = std::stoul( argv[3] );

  if( argc >= 5 )
    order = static_cast<unsigned>( std::stoul( argv[4] ) );

  // Data descriptor ---------------------------------------------------

  std::cerr << "* Calculating eccentricity data descriptor of order " << order << "...";

  auto eccentricity
    = eccentricities<Distance>( pointCloud, order );

  std::cerr << "finished\n";

  // Vietoris--Rips complex --------------------------------------------

  std::cerr << "* Calculating Vietoris--Rips complex with eps=" << epsilon << " and d=" << dimension << "...";

  FLANN<PointCloud, Distance> flannWrapper( pointCloud );

  auto K
    = buildVietorisRipsComplex( flannWrapper,
                                epsilon,
                                unsigned( dimension ),
                                eccentricity.begin(), eccentricity.end() );

  std::cerr << "finished\n"
            << "* Obtained simplicial complex with " << K.size() << " simplices\n";

  // Persistent homology -----------------------------------------------

  std::cerr << "* Calculating persistence diagrams...";

  auto diagrams
    = calculatePersistenceDiagrams( K );

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
