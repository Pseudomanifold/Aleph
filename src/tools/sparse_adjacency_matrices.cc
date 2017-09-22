#include <aleph/geometry/RipsExpander.hh>

#include <aleph/math/KahanSummation.hh>

#include <aleph/persistentHomology/Calculation.hh>

#include <aleph/topology/FloydWarshall.hh>
#include <aleph/topology/Simplex.hh>
#include <aleph/topology/SimplicialComplex.hh>

#include <aleph/topology/filtrations/Data.hh>
#include <aleph/topology/filtrations/Degree.hh>

#include <aleph/topology/io/GML.hh>
#include <aleph/topology/io/SparseAdjacencyMatrix.hh>

#include <aleph/utilities/Filesystem.hh>
#include <aleph/utilities/Format.hh>

#include <getopt.h>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <cmath>

using DataType          = float;
using VertexType        = std::size_t;
using Simplex           = aleph::topology::Simplex<DataType, VertexType>;
using SimplicialComplex = aleph::topology::SimplicialComplex<Simplex>;

std::vector<DataType> closenessCentrality( const SimplicialComplex& K )
{
  auto M = aleph::topology::floydWarshall( K, 1 );
  auto n = M.numRows();

  std::vector<DataType> result;
  result.reserve( n );

  for( decltype(n) i = 0; i < n; i++ )
  {
    aleph::math::KahanSummation<DataType> sum = DataType();

    for( decltype(n) j = 0; j < n; j++ )
      if( std::isfinite( M(i,j) ) )
        sum += M(i,j);

    result.push_back( DataType(n) / sum );
  }

  return result;
}

int main( int argc, char** argv )
{
  static option commandLineOptions[] =
  {
    { "closeness-centrality", no_argument, nullptr, 'c' },
    { "sum"                 , no_argument, nullptr, 's' },
    { nullptr               , 0          , nullptr,  0  }
  };

  bool calculateClosenessCentrality = false;
  bool useSumOfDegrees              = false;

  {
    int option = 0;
    while( ( option = getopt_long( argc, argv, "cs", commandLineOptions, nullptr ) ) != -1 )
    {
      switch( option )
      {
      case 'c':
        calculateClosenessCentrality = true;
        break;
      case 's':
        useSumOfDegrees = true;
        break;
      }
    }
  }

  if( ( argc - optind ) < 1 )
    return -1;

  std::string filename = argv[1];

  std::vector<SimplicialComplex> simplicialComplexes;
  std::vector<std::string> labels;

  aleph::topology::io::SparseAdjacencyMatrixReader reader;
  reader.setReadGraphLabels();

  std::cerr << "* Reading '" << filename << "'...";

  reader( filename, simplicialComplexes );

  std::cerr << "finished\n"
            << "* Read " << simplicialComplexes.size() << " simplicial complexes\n";

  // Calculate closeness centrality ------------------------------------

  if( calculateClosenessCentrality )
  {
    std::size_t index = 0;

    for( auto&& K : simplicialComplexes )
    {
      K.sort();
      auto cc     = closenessCentrality( K );
      auto output = "/tmp/"
                    + aleph::utilities::format( index, simplicialComplexes.size() )
                    + "_closeness_centrality.txt";

      std::ofstream out( output );
      for( auto&& value : cc )
        out << value << "\n";

      ++index;
    }
  }

  // Expand simplicial complexes ---------------------------------------

  aleph::geometry::RipsExpander<SimplicialComplex> expander;

  // TODO: make expansion configurable; does it make sense to expand the
  // complexes at all?

  // Calculate degrees -------------------------------------------------

  DataType maxDegree = 0;

  std::cerr << "* Calculating degree-based filtration...";

  for( auto&& K : simplicialComplexes )
  {
    std::vector<DataType> degrees;
    aleph::topology::filtrations::degrees( K, std::back_inserter( degrees ) );

    // TODO: check degree for isolated vertices?

    if( !degrees.empty() )
    {
      maxDegree
        = std::max( maxDegree,
                    *std::max_element( degrees.begin(), degrees.end() ) );
    }

    if( useSumOfDegrees )
      K = expander.assignData( K, degrees.begin(), degrees.end(), DataType(0), [] ( DataType a, DataType b ) { return a+b; } );
    else
      K = expander.assignMaximumData( K, degrees.begin(), degrees.end() );

    K.sort( aleph::topology::filtrations::Data<Simplex>() );
  }

  std::cerr << "finished\n"
            << "* Identified maximum degree as D=" << maxDegree << "\n";

  // Store graphs ------------------------------------------------------

  {
    aleph::topology::io::GMLWriter writer;

    for( std::size_t i = 0; i < simplicialComplexes.size(); i++ )
    {
      auto filename = "/tmp/"
                      + aleph::utilities::format( i, simplicialComplexes.size() )
                      + ".gml";

      std::cerr << "* Storing graph in '" << filename << "'...";

      writer( filename, simplicialComplexes[i] );

      std::cerr << "finished\n";
    }
  }

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
            out << point.x() << "\t" << 2 * maxDegree << "\n";
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
}
