#include <tests/Base.hh>

#include <aleph/persistentHomology/PhiPersistence.hh>

#include <aleph/topology/Simplex.hh>
#include <aleph/topology/SimplicialComplex.hh>

#include <vector>

template <class T> void test()
{
  using Simplex           = aleph::topology::Simplex<T>;
  using SimplicialComplex = aleph::topology::SimplicialComplex<Simplex>;

  std::vector<Simplex> simplices =
  {
    {0}, {1}, {2}, {3}, {4}
  };

  SimplicialComplex K( simplices.begin(), simplices.end() );

  // TODO: add real function
  aleph::partition(K, [] () {} );
}

int main(int, char**)
{
}
