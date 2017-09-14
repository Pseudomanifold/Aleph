#include <aleph/persistentHomology/ConnectedComponents.hh>

#include <aleph/topology/io/SimplicialComplexReader.hh>

#include <aleph/topology/Simplex.hh>
#include <aleph/topology/SimplicialComplex.hh>

#include <iostream>
#include <string>
#include <unordered_map>

using DataType          = double;
using VertexType        = unsigned;
using Simplex           = aleph::topology::Simplex<DataType, VertexType>;
using SimplicialComplex = aleph::topology::SimplicialComplex<Simplex>;

namespace
{

struct SizeFunctor
{
  void initialize( VertexType v )
  {
    _componentSize[v] = 1;
  }

  void operator()( VertexType younger,
                   VertexType older,
                   DataType creation,
                   DataType destruction )
  {
    _componentSize[older] += _componentSize[younger];

    std::cerr << "DEBUG:\n"
              << "  MERGE " << younger << " -> " << older << "\n"
              << "  PAIR  " << "(" << creation << "," << destruction << ")\n";
  }

  void operator()( VertexType root,
                   DataType creation )
  {
    std::cerr << "DEBUG:\n"
              << " ROOT     " << root << "\n"
              << " CREATION " << creation << "\n";
  }

  std::unordered_map<VertexType, unsigned> _componentSize;
};

}

int main( int argc, char** argv )
{
  if( argc <= 1 )
    return -1;

  std::string filename = argv[1];

  SimplicialComplex K;

  aleph::topology::io::SimplicialComplexReader reader;
  reader( filename, K );

  SizeFunctor sf;

  auto&& tuple = aleph::calculateZeroDimensionalPersistenceDiagram<Simplex, aleph::traits::PersistencePairingCalculation<aleph::PersistencePairing<VertexType> > >( K, sf );
  auto&& pd    = std::get<0>( tuple );
  auto&& pp    = std::get<1>( tuple );

  std::cout << pd << "\n";
}
