#include "tests/Base.hh"

#include "topology/Mesh.hh"
#include "topology/MorseSmaleComplex.hh"

#include <vector>

void test1()
{
  aleph::topology::Mesh<double> M;

  M.addVertex( 0.0, 0.0, 0.0 );
  M.addVertex( 0.0, 1.0, 0.0 );
  M.addVertex( 1.0, 0.0, 0.0 );
  M.addVertex( 1.5, 1.0, 0.0 );

  std::vector<unsigned> f1 = { 0, 1, 2 };
  std::vector<unsigned> f2 = { 2, 1, 3 };

  M.addFace( f1.begin(), f1.end() );
  M.addFace( f2.begin(), f2.end() );
}

void test2()
{
  aleph::topology::Mesh<double> M;

  M.addVertex( 0.0, 0.0, 0.0, 0.0 );
  M.addVertex( 1.0, 0.0, 0.0, 1.0 );
  M.addVertex( 2.0, 0.0, 0.0, 0.0 );
  M.addVertex( 0.0, 1.0, 0.0, 1.0 );
  M.addVertex( 1.0, 1.0, 0.0, 2.0 );
  M.addVertex( 2.0, 1.0, 0.0, 1.0 );
  M.addVertex( 0.0, 2.0, 0.0, 0.0 );
  M.addVertex( 1.0, 2.0, 0.0, 1.0 );
  M.addVertex( 2.0, 2.0, 0.0, 0.0 );

  std::vector<unsigned> f1 = { 0, 1, 4, 3 };
  std::vector<unsigned> f2 = { 1, 2, 5, 4 };
  std::vector<unsigned> f3 = { 4, 5, 8, 7 };
  std::vector<unsigned> f4 = { 3, 4, 7, 6 };

  M.addFace( f1.begin(), f1.end() );
  M.addFace( f2.begin(), f2.end() );
  M.addFace( f3.begin(), f3.end() );
  M.addFace( f4.begin(), f4.end() );

  aleph::topology::MorseSmaleComplex<decltype(M)> msc;
  msc( M );
}

int main(int, char**)
{
  test1();
  test2();
}
