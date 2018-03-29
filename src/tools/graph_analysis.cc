/**
  @file  graph_analysis.cc
  @brief Tool for calculating persistent homology of graphs

  This is a tool shipped by 'Aleph - A Library for Exploring Persistent
  Homology'. It uses a simple degree filtration to convert a graph that
  is specified as a list of edges into a simplicial complex.

  Original author: Bastian Rieck
*/

#include <aleph/geometry/RipsExpander.hh>

#include <aleph/persistentHomology/Calculation.hh>

#include <aleph/topology/Simplex.hh>
#include <aleph/topology/SimplicialComplex.hh>

#include <aleph/topology/filtrations/Data.hh>
#include <aleph/topology/filtrations/Degree.hh>

#include <aleph/topology/io/EdgeLists.hh>

#include <aleph/utilities/String.hh>

#include <getopt.h>

#include <iostream>
#include <vector>

using DataType          = unsigned;
using VertexType        = unsigned;
using Simplex           = aleph::topology::Simplex<DataType, VertexType>;
using SimplicialComplex = aleph::topology::SimplicialComplex<Simplex>;

int main( int argc, char** argv )
{
  static option commandLineOptions[] =
  {
    { "infinity"            , required_argument, nullptr, 'f' },
    { "loops"               , no_argument      , nullptr, 'l' },
    { nullptr               , 0                , nullptr,  0  }
  };

  DataType infinity   = DataType();
  bool calculateLoops = false;

  {
    int option = 0;
    while( ( option = getopt_long( argc, argv, "f:l", commandLineOptions, nullptr ) ) != -1 )
    {
      switch( option )
      {
      case 'f':
        infinity = aleph::utilities::convert<DataType>( optarg );
        break;
      case 'l':
        calculateLoops = true;
        break;
      }
    }
  }

  if( ( argc - optind ) < 1 )
    return -1;

  std::string filename = argv[optind++];

  SimplicialComplex K;

  aleph::topology::io::EdgeListReader reader;
  reader.setReadWeights();
  reader.setTrimLines();

  std::cerr << "* Reading '" << filename << "'...";

  reader( filename, K );

  std::cerr << "finished\n"
            << "* Read simplicial complex with " << K.size() << " simplices\n";


  // Calculate degrees -------------------------------------------------

  std::cerr << "* Calculating degree-based filtration...";

  aleph::geometry::RipsExpander<SimplicialComplex> expander;

  {
    std::vector<unsigned> degrees;
    aleph::topology::filtrations::degrees( K, std::back_inserter( degrees ) );

    K = expander.assignMaximumData( K, degrees.begin(), degrees.end() );
    K.sort( aleph::topology::filtrations::Data<Simplex>() );
  }

  std::cerr << "finished\n";

#if 0
  // Calculate persistent homology -------------------------------------

  {
    std::size_t index = 0;

    for( auto&& K : simplicialComplexes )
    {
      bool dualize                    = true;
      bool includeAllUnpairedCreators = true;

      auto diagrams
        = aleph::calculatePersistenceDiagrams( K,
                                               dualize,
                                               includeAllUnpairedCreators );

      for( auto&& diagram : diagrams )
      {
        diagram.removeDiagonal();

        auto output = "/tmp/"
                      + aleph::utilities::format( index, simplicialComplexes.size() )
                      + "_d"
                      + std::to_string( diagram.dimension() )
                      + ".txt";

        std::ofstream out( output );

        for( auto&& point : diagram )
        {
          if( point.isUnpaired() )
            out << point.x() << "\t" << infinity * maxDegree << "\n";
          else
            out << point.x() << "\t" << point.y() << "\n";
        }
      }

      ++index;
    }
  }

  // Store labels ------------------------------------------------------

  {
    std::vector<std::string> labels;
    reader.graphLabels( std::back_inserter( labels ) );

    std::ofstream out( "/tmp/"
                      + aleph::utilities::basename( filename )
                      + ".txt" );

    for( auto&& label : labels )
      out << label << "\n";
  }
#endif
}
