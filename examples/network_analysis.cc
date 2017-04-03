#include "geometry/RipsExpander.hh"

#include "filtrations/Data.hh"

#include "topology/io/EdgeLists.hh"

#include "persistentHomology/Calculation.hh"

#include "topology/Simplex.hh"
#include "topology/SimplicialComplex.hh"

#include <iostream>
#include <fstream>
#include <string>

using DataType          = double;
using VertexType        = unsigned;

using Simplex           = aleph::topology::Simplex<DataType, VertexType>;
using SimplicialComplex = aleph::topology::SimplicialComplex<Simplex>;

int main( int argc, char** argv )
{
  if( argc <= 1 )
    return -1;

  aleph::topology::io::EdgeListReader reader;

  std::string input = argv[1];
  std::ifstream in( input );

  std::cerr << "* Reading '" << input << "'...";

  SimplicialComplex K;
  reader( in, K );

  std::cerr << "finished\n"
            << "* Extracted simplicial complex with " << K.size() << " simplices\n";

  std::cerr << "* Expanding simplicial complex...";

  // TODO: Make expansion configurable
  aleph::geometry::RipsExpander<SimplicialComplex> ripsExpander;
  K = ripsExpander( K, 2 );
  K = ripsExpander.assignMaximumWeight( K );

  std::cerr << "...finished\n"
            << "* Expanded complex has " << K.size() << " simplices\n";

  std::cerr << "* Establishing filtration order...";

  K.sort(
    aleph::filtrations::Data<Simplex>()
  );

  std::cerr << "...finished\n";

  std::cerr << "* Calculating persistent homology...";

  auto diagrams
    = aleph::calculatePersistenceDiagrams( K );

  std::cerr << "...finished\n";

  for( auto&& D : diagrams )
  {
    D.removeDiagonal();

    std::cout << "# Persistence diagram <" << input << ">\n"
              << "#\n"
              << "# Dimension   : " << D.dimension() << "\n"
              << "# Entries     : " << D.size() << "\n"
              << "# Betti number: " << D.betti() << "\n"
              << D << "\n\n";
  }
}
