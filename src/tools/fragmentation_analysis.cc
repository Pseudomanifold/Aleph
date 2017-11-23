#include <aleph/topology/io/GML.hh>

#include <aleph/persistentHomology/Calculation.hh>

#include <aleph/topology/Simplex.hh>
#include <aleph/topology/SimplicialComplex.hh>

#include <aleph/topology/filtrations/Data.hh>

#include <iostream>
#include <sstream>
#include <string>

int main( int argc, char** argv )
{
  if( argc <= 1 )
    return -1;

  std::string filename = argv[1];

  using DataType          = double;
  using VertexType        = unsigned short;
  using Simplex           = aleph::topology::Simplex<DataType, VertexType>;
  using SimplicialComplex = aleph::topology::SimplicialComplex<Simplex>;

  // Input -------------------------------------------------------------

  SimplicialComplex K;

  std::cerr << "* Reading " << filename << "...";

  {
    aleph::topology::io::GMLReader reader;
    reader( filename, K );
  }

  std::cerr << "finished\n";

  // Data assignment ---------------------------------------------------
  //
  // TODO: assign selected data attribute as the weight of the
  // simplicial complex, either for vertices or for edges

  // Filtration --------------------------------------------------------

  K.sort(
    aleph::topology::filtrations::Data<Simplex>()
  );

  // Persistent homlogy ------------------------------------------------

  std::cerr << "* Calculating persistent homology...";

  bool dualize                    = true;
  bool includeAllUnpairedCreators = true;

  auto diagrams
    = aleph::calculatePersistenceDiagrams( K,
                                           dualize,
                                           includeAllUnpairedCreators );

  std::cerr << "finished\n";

  for( auto&& D : diagrams )
  {
    D.removeDiagonal();

    std::ostringstream stream;

    stream << "# Persistence diagram <" << filename << ">\n"
           << "#\n"
           << "# Dimension   : " << D.dimension() << "\n"
           << "# Entries     : " << D.size() << "\n"
           << "# Betti number: " << D.betti() << "\n"
           << D;

    if( D != diagrams.back() )
      stream << "\n\n";

    std::cout << stream.str();
  }
}
