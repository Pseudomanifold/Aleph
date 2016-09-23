#include "distances/detail/Matrix.hh"
#include "distances/detail/Munkres.hh"

using namespace aleph::distances::detail;

int main()
{
  Matrix<unsigned short> m( 3 );

  m(0,0) = 1; m(0,1) = 2; m(0,2) = 3;
  m(1,0) = 2; m(1,1) = 4; m(1,2) = 6;
  m(2,0) = 3; m(2,1) = 6; m(2,2) = 9;

  Munkres<unsigned short> solver( m );
  solver();
}
