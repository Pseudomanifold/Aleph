/**
  This is a tool shipped by 'Aleph - A Library for Exploring Persistent
  Homology'.

  Its purpose is to calculate zero-dimensional persistence diagrams of
  spectra. This is supposed to yield a simple feature descriptor which
  in turn might be used in machine learning pipepline.

  Input:  filename
  Output: persistence diagram

  The persistence diagram represents the superlevel set filtration of
  the input data. This permits us to quantify the number of maxima in
  a data set.

  Original author: Bastian Rieck
*/

#include <aleph/persistenceDiagrams/PersistenceDiagram.hh>

#include <aleph/persistentHomology/Calculation.hh>

#include <aleph/topology/io/FlexSpectrum.hh>

#include <aleph/topology/Simplex.hh>
#include <aleph/topology/SimplicialComplex.hh>

#include <algorithm>
#include <iostream>
#include <string>

#include <cassert>

using DataType           = unsigned;
using VertexType         = unsigned;
using Simplex            = aleph::topology::Simplex<DataType, VertexType>;
using SimplicialComplex  = aleph::topology::SimplicialComplex<Simplex>;
using PersistenceDiagram = aleph::PersistenceDiagram<DataType>;

int main( int argc, char** argv )
{
  if( argc <= 1 )
    return -1;

  std::string input = argv[1];

  // Parse input -------------------------------------------------------

  std::cerr << "* Reading '" << input << "'...";

  SimplicialComplex K;

  aleph::topology::io::FlexSpectrumReader reader;
  reader( input, K );

  std::cerr << "finished\n";

  // Calculate persistent homology -------------------------------------

  std::cerr << "* Calculating persistent homology...";

  auto diagrams = aleph::calculatePersistenceDiagrams( K );

  // Need to ensure that we are actually doing the right thing here, so
  // I rather check *everything* that might go awry.
  assert( diagrams.empty() == false );
  assert( diagrams.size()  == 1 );

  std::cerr << "finished\n";

  // Output ------------------------------------------------------------

  if( diagrams.empty() == false )
  {
    auto&& D = diagrams.front();

    assert( D.dimension() == 0 );
    assert( D.betti()     == 1 );

    D.removeDiagonal();

    // This ensures that the global maximum is paired with the global
    // minimum of the persistence diagram. This is valid because each
    // function has finite support and is bounded from below.
    std::transform( D.begin(), D.end(), D.begin(),
      [] ( const PersistenceDiagram::Point& p )
      {
        // TODO: we should check whether zero is really the smallest
        // value
        if( p.isUnpaired() )
          return PersistenceDiagram::Point( p.x(), DataType() );
        else
          return PersistenceDiagram::Point( p );
      }
    );

    std::cout << D << "\n";
  }
}
