#include <aleph/persistentHomology/ConnectedComponents.hh>

#include <aleph/topology/filtrations/Data.hh>

#include <aleph/topology/io/SimplicialComplexReader.hh>

#include <aleph/topology/Simplex.hh>
#include <aleph/topology/SimplicialComplex.hh>

#include <iostream>
#include <map>
#include <string>
#include <utility>
#include <unordered_map>

using DataType          = double;
using VertexType        = unsigned;
using Simplex           = aleph::topology::Simplex<DataType, VertexType>;
using SimplicialComplex = aleph::topology::SimplicialComplex<Simplex>;
using Filtration        = aleph::topology::filtrations::Data<Simplex>;
using Edge              = std::pair<VertexType, VertexType>;

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
                   DataType destruction,
                   VertexType u,
                   VertexType v )
  {
    if( u > v )
      std::swap(u,v);

    _componentSize[older]                 += _componentSize[younger];
    _edgeRelevance[ std::make_pair(u,v) ]  = _componentSize[older];

    std::cerr << "* Edge " << "(" << u << "," << v << "): "
              << _edgeRelevance[ std::make_pair(u,v) ]
              << " "
              << "[" << creation << "," << destruction << "]"
              << "\n";
  }

  void operator()( VertexType /* root */,
                   DataType /* creation */ )
  {
  }

  std::unordered_map<VertexType, unsigned> _componentSize;
  std::map<Edge, unsigned>                 _edgeRelevance;
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

  K.sort( Filtration() );

  SizeFunctor sf;

  auto&& tuple = aleph::calculateZeroDimensionalPersistenceDiagram<Simplex, aleph::traits::PersistencePairingCalculation<aleph::PersistencePairing<VertexType> > >( K, sf );
  auto&& pd    = std::get<0>( tuple );
  auto&& pp    = std::get<1>( tuple );

  std::cout << pd << "\n";
}
