#include "complexes/RipsExpander.hh"

#include "filtrations/Data.hh"

#include "io/EdgeLists.hh"

#include <iostream>
#include <fstream>
#include <string>

#include "Simplex.hh"
#include "SimplicialComplex.hh"

#include "PersistentHomologyCalculation.hh"

using DataType          = double;
using VertexType        = unsigned;

using Simplex           = aleph::Simplex<DataType, VertexType>;
using SimplicialComplex = aleph::SimplicialComplex<Simplex>;

int main( int argc, char** argv )
{
  if( argc <= 1 )
    return -1;

  aleph::io::EdgeListReader reader;

  std::string input = argv[1];
  std::ifstream in( input );

  std::cerr << "* Reading '" << input << "'...";

  SimplicialComplex K;
  reader( in, K );

  std::cerr << "finished\n"
            << "* Extracted simplicial complex with " << K.size() << " simplices\n";

  std::cerr << "* Expanding simplicial complex...";

  // TODO: Make expansion configurable
  aleph::complexes::RipsExpander<Simplex> ripsExpander;
  K = ripsExpander( K, 2 );

  std::cerr << "...finished\n"
            << "* Expanded complex has " << K.size() << " simplices\n";

  std::cerr << "* Establishing filtration order...";

  K.sort(
    aleph::filtrations::Data<Simplex>()
  );

  std::cerr << "...finished\n";

  std::cerr << "* Calculating persistent homology...";

  auto diagrams
    = calculatePersistenceDiagrams( K );

  std::cerr << "...finished\n";

  for( auto&& D : diagrams )
    std::cout << D << "\n";
}
