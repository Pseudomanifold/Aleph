#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <sstream>

#include <cmath>

#include "filtrations/Data.hh"

#include "geometry/RipsExpander.hh"

#include "persistentHomology/ConnectedComponents.hh"

#include "topology/CliqueGraph.hh"
#include "topology/ConnectedComponents.hh"
#include "topology/Simplex.hh"
#include "topology/SimplicialComplex.hh"

#include "topology/io/EdgeLists.hh"
#include "topology/io/GML.hh"

#include "utilities/Filesystem.hh"

using DataType           = double;
using VertexType         = unsigned;
using Simplex            = aleph::topology::Simplex<DataType, VertexType>;
using SimplicialComplex  = aleph::topology::SimplicialComplex<Simplex>;
using PersistenceDiagram = aleph::PersistenceDiagram<DataType>;

SimplicialComplex filterSimplicialComplex( const SimplicialComplex& K, DataType threshold )
{
  std::vector<Simplex> simplices;

  std::copy_if( K.begin(), K.end(), std::back_inserter( simplices ),
                [&threshold] ( const Simplex& s )
                {
                  // TODO: Less than or equal vs. less than?
                  return s.data() < threshold;
                } );

  return SimplicialComplex( simplices.begin(), simplices.end() );
}

std::string formatOutput( const std::string& prefix, unsigned k, unsigned K )
{
  std::ostringstream stream;
  stream << prefix;
  stream << std::setw( int( std::log10( K ) + 1 ) ) << std::setfill( '0' ) << k;
  stream << ".txt";

  return stream.str();
}

void usage()
{
  std::cerr << "Usage: clique-persistence-diagram FILE K\n"
            << "\n"
            << "Calculates the clique persistence diagram for FILE, which is\n"
            << "supposed to be a weighted graph. The K parameter denotes the\n"
            << "maximum dimension of a simplex for extracting a clique graph\n"
            << "and tracking persistence of clique communities.\n\n";
}

int main( int argc, char** argv )
{
  if( argc <= 2 )
  {
    usage();
    return -1;
  }

  std::string filename = argv[1];
  unsigned maxK        = static_cast<unsigned>( std::stoul( argv[2] ) );

  SimplicialComplex K;

  // Input -------------------------------------------------------------

  std::cerr << "* Reading '" << filename << "'...";

  if( aleph::utilities::extension( filename ) == ".gml" )
  {
    aleph::topology::io::GMLReader reader;
    reader( filename, K );
  }
  else
  {
    aleph::io::EdgeListReader reader;
    reader.setReadWeights( true );
    reader.setTrimLines( true );

    reader( filename, K );
  }

  std::cerr << "finished\n";

  DataType maxWeight = std::numeric_limits<DataType>::lowest();
  for( auto&& simplex : K )
    maxWeight = std::max( maxWeight, simplex.data() );

  // Expansion ---------------------------------------------------------

  aleph::geometry::RipsExpander<SimplicialComplex> ripsExpander;
  K = ripsExpander( K, maxK );
  K = ripsExpander.assignMaximumWeight( K );

  K.sort( aleph::filtrations::Data<Simplex>() );

  // Stores the accumulated persistence of vertices. Persistence
  // accumulates if a vertex participates in a clique community.
  std::map<VertexType, double> accumulatedPersistenceMap;

  for( unsigned k = 1; k <= maxK; k++ )
  {
    std::cerr << "* Extracting " << k << "-cliques graph...";

    auto C
        = aleph::topology::getCliqueGraph( K, k );

    C.sort( aleph::filtrations::Data<Simplex>() );

    std::cerr << "finished\n";

    std::cerr << "* " << k << "-cliques graph has " << C.size() << " simplices\n";

    if( C.empty() )
    {
      std::cerr << "* Stopping here because no further cliques for processing exist\n";
      break;
    }

    auto pd
        = aleph::calculateZeroDimensionalPersistenceDiagram( C );

    // TODO: What about duplicates?
    for( auto&& point : pd )
    {
      if( point.x() == point.y() )
        continue;

      // TODO:
      //   - Extract simplicial complex of given threshold
      //   - Calculate connected components
      //   - Extract connected component containing the creator simplex
      //   - Traverse the connected component and assign persistence

      auto epsilon         = point.y();
      auto filteredComplex = filterSimplicialComplex( C, epsilon );
      auto uf              = calculateConnectedComponents( filteredComplex );

      std::set<VertexType> roots;
      uf.roots( std::inserter( roots, roots.begin() ) );

      std::set<VertexType> cliqueVertices;

      for( auto&& root : roots )
      {
        // Only consider roots that 'fit' the current creation threshold
        // value. In the filtered complex, other connected components of
        // a different persistence may still exist.
        //
        // TODO: Does it make sense to filter with upper _and_ lower
        // bounds?
        if( C.find( root )->data() == point.x() )
        {
          std::vector<VertexType> vertices;
          uf.get( root, std::back_inserter( vertices ) );

          for( auto&& vertex : vertices )
          {
            auto it = filteredComplex.find( vertex );
            if( it == filteredComplex.end() )
              throw std::runtime_error( "Unable to identify simplex in filtered simplicial complex" );

            cliqueVertices.insert( it->begin(), it->end() );
          }
        }
      }

      for( auto&& cliqueVertex : cliqueVertices )
        accumulatedPersistenceMap[cliqueVertex] += std::isfinite( point.persistence() ) ? point.persistence() : maxWeight - point.x();
    }

    {
      using namespace aleph::utilities;
      auto outputFilename = formatOutput( "/tmp/" + stem( basename( filename ) ) + "_k", k, maxK );

      std::cerr << "* Storing output in '" << outputFilename << "'...\n";

      pd.removeDiagonal();

      std::transform( pd.begin(), pd.end(), pd.begin(),
                      [&maxWeight] ( const PersistenceDiagram::Point& p )
                      {
                        if( !std::isfinite( p.y() ) )
                          return PersistenceDiagram::Point( p.x(), maxWeight );
                        else
                          return p;
                      } );

      std::ofstream out( outputFilename );
      out << "# Original filename: " << filename << "\n";
      out << "# k                : " << k        << "\n";
      out << pd << "\n";
    }
  }

  {
    using namespace aleph::utilities;
    auto outputFilename = "/tmp/" + stem( basename( filename ) ) + ".txt";

    std::cerr << "* Storing accumulated persistence values in '" << outputFilename << "'...\n";

    std::ofstream out( outputFilename );

    for( auto&& pair : accumulatedPersistenceMap )
      out << pair.first << "\t" << pair.second << "\n";
  }
}
