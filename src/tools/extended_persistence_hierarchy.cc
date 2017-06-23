/*
  This is a tool shipped by 'Aleph - A Library for Exploring Persistent
  Homology'.

  It loads VTK files (structured grids) or TXT (1D functions) files and
  calculates the extended persistence hierarchy described in:

    Hierarchies and Ranks for Persistence Pairs
    Bastian Rieck, Heike Leitte, and Filip Sadlo
    Proceedings of TopoInVis 2017, Japan

  The output of the tool is a list of nodes for the hierarchy, followed
  by a list of edges. Each node is identified by an ID, followed by its
  corresponding persistence pair entry. An edge consists of two node ID
  values, connected via "--".

  The following snippet gives a small demonstration:

    0: 0 infty
    1: 1 2
    2: 3 4

    0 -- 1
    0 -- 2

  This output may subsequently be used in an auxiliary Python script in
  order to analyse and compare different hierarchies.

  TODO: Document this
*/

#include "persistenceDiagrams/Calculation.hh"
#include "persistenceDiagrams/PersistenceDiagram.hh"

#include "persistentHomology/ExtendedPersistenceHierarchy.hh"
#include "persistentHomology/PersistencePairing.hh"

#include "topology/Simplex.hh"
#include "topology/SimplicialComplex.hh"
#include "topology/UnionFind.hh"

#include "topology/filtrations/Data.hh"

#include "topology/io/Function.hh"
#include "topology/io/VTK.hh"

#include "utilities/Filesystem.hh"

#include <algorithm>
#include <functional>
#include <fstream>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

// TODO: Replace this as soon as possible with a more modern option
// parser interface.
#include <getopt.h>

using DataType           = double;
using VertexType         = unsigned;
using Simplex            = aleph::topology::Simplex<DataType, VertexType>;
using SimplicialComplex  = aleph::topology::SimplicialComplex<Simplex>;
using PersistenceDiagram = aleph::PersistenceDiagram<DataType>;

void usage()
{
  std::cerr << "Usage: extended_persistence_hierarchy [--superlevels] [--sublevels] FILES\n"
            << "\n"
            << "Calculates the extended persistence hierarchy of a set of VTK files or 1D\n"
            << "functions stored in FILES. By default, a filtration based on the sublevel\n"
            << "sets is used. This may either be enforced or modified by using one of the\n"
            << "long options specified above.\n"
            << "\n"
            << "The hierarchy is written to STDOUT.\n"
            << "\n"
            << "Flags:\n"
            << "  -s: use sublevel set filtration\n"
            << "  -S: use superlevel set filtration\n"
            << "\n";
}

int main( int argc, char** argv )
{
  static option commandLineOptions[] =
  {
    { "superlevels", no_argument, nullptr, 'S' },
    { "sublevels"  , no_argument, nullptr, 's' },
    { nullptr      , 0          , nullptr,  0  }
  };

  bool calculateSuperlevelSets = false;

  int option = 0;
  while( ( option = getopt_long( argc, argv, "Ss", commandLineOptions, nullptr ) ) != -1 )
  {
    switch( option )
    {
    case 'S':
      calculateSuperlevelSets = true;
      break;
    case 's':
      calculateSuperlevelSets = false;
      break;
    default:
      break;
    }
  }

  if( ( argc - optind ) < 1 )
  {
    usage();
    return -1;
  }

  std::vector<std::string> filenames;
  filenames.reserve( argc - optind );

  for( int i = optind; i < argc; i++ )
    filenames.push_back( argv[i] );

  std::vector<SimplicialComplex> simplicialComplexes;
  simplicialComplexes.reserve( filenames.size() );

  // Select a functor for calculating the weights when reading
  // a simplicial complex.
  auto functor = calculateSuperlevelSets
               ? [] ( DataType a, DataType b ) { return std::min(a,b); }
               : [] ( DataType a, DataType b ) { return std::max(a,b); };

  for( auto&& filename : filenames )
  {
    std::cerr << "* Reading '" << filename << "'...";

    if( aleph::utilities::extension( filename ) == ".vtk" )
    {
      SimplicialComplex K;

      aleph::topology::io::VTKStructuredGridReader reader;
      reader( filename, K, functor );

      if( calculateSuperlevelSets )
        K.sort( aleph::topology::filtrations::Data<Simplex, std::greater<DataType> >() );
      else
        K.sort( aleph::topology::filtrations::Data<Simplex, std::less<DataType> >() );

      simplicialComplexes.emplace_back( K );
    }
    else
    {
      auto complexes
        = aleph::topology::io::loadFunctions<SimplicialComplex>( filename,
                                                                 functor );

      for( auto&& K : complexes )
      {
        if( calculateSuperlevelSets )
          K.sort( aleph::topology::filtrations::Data<Simplex, std::greater<DataType> >() );
        else
          K.sort( aleph::topology::filtrations::Data<Simplex, std::less<DataType> >() );
      }

      simplicialComplexes.insert( simplicialComplexes.end(),
                                  complexes.begin(), complexes.end() );
    }

    std::cerr << "finished\n";
  }

  std::vector<PersistenceDiagram> persistenceDiagrams;
  persistenceDiagrams.reserve( simplicialComplexes.size() );

  for( auto&& K : simplicialComplexes )
  {
    aleph::ExtendedPersistenceHierarchy<Simplex> eph;
    auto ppe                = eph( K );   // persistence pairing & edges
    auto persistencePairing = ppe.first;  // persistence pairing
    auto edges              = ppe.second; // edges

    // Calculate the corresponding persistence diagram, clean it
    // perfunctorily, and finally sort it by increasing creation
    // threshold. Subsequently, we will check which diagrams are
    // unique.
    {
      auto diagrams = aleph::makePersistenceDiagrams( persistencePairing, K );
      auto diagram  = diagrams.front();

      diagram.removeDiagonal();

      std::sort( diagram.begin(), diagram.end(), [] ( const typename PersistenceDiagram::Point& a, const typename PersistenceDiagram::Point& b )
                                                 {
                                                   return a.x() < b.x() || ( a.x() == b.x() && a.y() < b.y() );
                                                 } );

      persistenceDiagrams.push_back( diagram );
    }

    // Enumerate all vertices in the hierarchy -------------------------

    std::set<VertexType> vertices;

    for( auto edge : edges )
    {
      vertices.insert( edge.first );
      vertices.insert( edge.second );
    }

    // Display nodes of the hierarchy ----------------------------------

    {
      unsigned index = 0;
      for( auto&& vertex : vertices )
      {
        auto itCreator = K.find( Simplex(vertex) );
        if( itCreator != K.end() )
        {
          auto itPair         = persistencePairing.find( static_cast<VertexType>( K.index( *itCreator ) ) );
          auto destroyerIndex = itPair->second;

          if( destroyerIndex < K.size() )
            std::cout << index << ": " << itCreator->data() << "\t" << K.at( destroyerIndex ).data() << "\n";
          else
            std::cout << index << ": " << itCreator->data() << "\t" << std::numeric_limits<DataType>::infinity() << "\n";
        }

        ++index;
      }

      std::cout << "\n";
    }

    // Display edges of the hierarchy ----------------------------------

    for( auto&& edge : edges )
    {
      auto uIndex = std::distance( vertices.begin(), vertices.find( edge.first ) );
      auto vIndex = std::distance( vertices.begin(), vertices.find( edge.second ) );

      std::cout << uIndex << " -- " << vIndex << "\n";
    }

    std::cout << "\n\n";
  }

  // Determine which persistence diagrams are unique among the
  // discovered ones...
  {
    std::vector<std::size_t> indices( persistenceDiagrams.size() );
    std::iota( indices.begin(), indices.end(), 0 );

    aleph::topology::UnionFind<std::size_t> uf( indices.begin(), indices.end() );

    for( std::size_t i = 0; i < persistenceDiagrams.size(); i++ )
    {
      auto&& D1 = persistenceDiagrams.at(i);

      for( std::size_t j = i+1; j < persistenceDiagrams.size(); j++ )
      {
        auto&& D2 = persistenceDiagrams.at(j);

        if( D1 == D2 )
          uf.merge( j, i );
      }
    }

    std::vector<std::size_t> roots;
    uf.roots( std::back_inserter( roots ) );

    std::ofstream out( "/tmp/classes.txt" );

    std::cerr << "* Out of " << persistenceDiagrams.size() << " persistence diagrams, there are " << roots.size() << " unique ones\n";

    unsigned C = 0;
    for( auto&& root : roots )
    {
      std::vector<std::size_t> indices;
      uf.get( root, std::back_inserter( indices ) );

      if( indices.size() > 1 )
      {
        std::cerr << "  Class " << C << ": " << indices.size() << "\n";

        C++; // ha ha

        for( auto&& index : indices )
          out << index << " ";

        out << "\n";
      }
    }
  }
}
