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

    attributes = reader.getEdgeAttributeNames();

    std::cerr << "* Available edge attributes: ";

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

  if( not edgeAttribute.empty() )
  {
    std::cerr << "* Using edge attribute '" << edgeAttribute << "' to assign weights...";

    // TODO: simplex traversal could be optimized by only taking edges
    // into account
    for( auto it = K.begin(); it != K.end(); ++it )
    {
      if( it->dimension() == 1 )
      {
        // Make a copy of the simplex because the container
        // configuration does not permit direct replacement
        // of any simplex.
        auto simplex = *it;
        auto source  = std::to_string( simplex[0] );
        auto target  = std::to_string( simplex[1] );

        // TODO: this is stupid and wasteful; the lookup could be
        // improved for larger graphs
        auto value1 = reader.getEdgeAttribute( source, target, edgeAttribute );
        auto value2 = reader.getEdgeAttribute( target, source, edgeAttribute );

        DataType w   = DataType();
        bool success = false;

        if( not value1.empty() )
          w = aleph::utilities::convert<DataType>( value1, success );
        else if( not value2.empty() )
          w = aleph::utilities::convert<DataType>( value2, success );
        else
          throw std::runtime_error( "Unable to find edge attribute" );

        if( not success )
          throw std::runtime_error( "Unable to convert edge attribute to data type" );

        simplex.setData( w );
        success = K.replace( it, simplex );

        if( !success )
          throw std::runtime_error( "Unable to replace simplex in simplicial complex" );
      }
    }

    // Recalculate all weights in the simplicial complex. This should
    // *not* be necessary for low-dimensional complexes, i.e. graphs,
    // but it might be useful later and I do not want this as a bug.
    bool useMaximum                  = true;
    bool skipOneDimensionalSimplices = true;

    K.recalculateWeights( useMaximum, skipOneDimensionalSimplices );

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
