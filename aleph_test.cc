#include "BoundaryMatrix.hh"
#include "Vector.hh"
#include "IO.hh"

#include <iostream>

using namespace aleph;
using namespace representations;

using I  = unsigned;
using V  = Vector<I>;
using BM = BoundaryMatrix<V>;

int main()
{
  auto M = load<BM>( "Triangle.txt" );

  std::cout << M << std::endl;
}
