#include <aleph/geometry/RipsExpander.hh>

#include <aleph/topology/io/GML.hh>

#include <aleph/persistentHomology/Calculation.hh>

#include <aleph/topology/Simplex.hh>
#include <aleph/topology/SimplicialComplex.hh>

#include <aleph/topology/filtrations/Data.hh>

#include <aleph/utilities/String.hh>

#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

#include <getopt.h>

int main( int argc, char** argv )
{
  if( argc <= 1 )
    return -1;

  std::string nodeAttribute;
  std::string edgeAttribute;

  {
    static option commandLineOptions[] =
    {
      { "node-attribute", required_argument, nullptr, 'n' },
      { "edge-attribute", required_argument, nullptr, 'e' },
      { nullptr         , 0                , nullptr,  0  }
    };

    int option = 0;
    while( ( option = getopt_long( argc, argv, "e:n:", commandLineOptions, nullptr ) ) != -1 )
    {
      switch( option )
      {
      case 'n':
        nodeAttribute = optarg;
        break;

      case 'e':
        edgeAttribute = optarg;
        break;

      default:
        break;
      }
    }
  }

  std::string filename = argv[ optind++ ];

  using DataType          = double;
  using VertexType        = unsigned short;
  using Simplex           = aleph::topology::Simplex<DataType, VertexType>;
  using SimplicialComplex = aleph::topology::SimplicialComplex<Simplex>;

  // Input -------------------------------------------------------------

  SimplicialComplex K;

  std::cerr << "* Reading " << filename << "...";

  aleph::topology::io::GMLReader reader;
  reader( filename, K );

  std::cerr << "finished\n";

  {
    auto attributes = reader.getNodeAttributeNames();

    std::cerr << "* Available node attributes: ";

    for( auto&& attribute : attributes )
    {
      std::cerr << attribute;
      if( attribute != attributes.back() )
        std::cerr << ", ";
    }

    std::cerr << "\n";
  }

  // Data assignment ---------------------------------------------------
  //
  // TODO: assign selected data attribute as the weight of the
  // simplicial complex, either for vertices or for edges

  if( not nodeAttribute.empty() )
  {
    std::cerr << "* Using node attribute '" << nodeAttribute << "' to assign weights...";

    auto id_to_index  = reader.id_to_index<VertexType>();
    auto attributeMap = reader.getNodeAttribute( nodeAttribute );

    // Assumption: attribute map contains exactly the same number of
    // entries than the number of vertices in the graph
    std::vector<DataType> data( attributeMap.size() );

    for( auto&& pair : attributeMap )
    {
      bool success = false;
      auto index   = id_to_index[pair.first];
      auto value   = aleph::utilities::convert<DataType>( pair.second, success );

      if( !success )
        throw std::runtime_error( "Unable to convert attribute value to data type" );

      data.at( index ) = value;
    }

    aleph::geometry::RipsExpander<SimplicialComplex> expander;
    K = expander.assignMaximumData( K, data.begin(), data.end() );

    std::cerr << "finished\n";
  }

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
