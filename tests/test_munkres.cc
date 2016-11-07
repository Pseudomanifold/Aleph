#include <iostream>
#include <string>

#include "distances/detail/Matrix.hh"
#include "distances/detail/Munkres.hh"

#include "tests/Base.hh"

using namespace aleph::distances::detail;

template <class T> void threeByThree()
{
  ALEPH_TEST_BEGIN( "Solving a three-by-three matrix" );

  Matrix<T> m( 3 );

  m(0,0) = static_cast<T>(1); m(0,1) = static_cast<T>(2); m(0,2) = static_cast<T>(3);
  m(1,0) = static_cast<T>(2); m(1,1) = static_cast<T>(4); m(1,2) = static_cast<T>(6);
  m(2,0) = static_cast<T>(3); m(2,1) = static_cast<T>(6); m(2,2) = static_cast<T>(9);

  Munkres<T> solver( m );
  solver();

  auto cost = solver.cost( m );

  ALEPH_ASSERT_THROW( cost  > T()   );
  ALEPH_ASSERT_THROW( cost == T(10) );

  ALEPH_TEST_END();
}

int main()
{
  threeByThree<short>         ();
  threeByThree<unsigned short>();
  threeByThree<int>           ();
  threeByThree<unsigned int>  ();
  threeByThree<long>          ();
  threeByThree<unsigned long> ();
  threeByThree<float>         ();
  threeByThree<double>        ();

#if 0

  {
    Matrix<unsigned short> m( 4 );

    m(0,0) = 82; m(0,1) = 83; m(0,2) = 69; m(0,3) = 92;
    m(1,0) = 77; m(1,1) = 37; m(1,2) = 49; m(1,3) = 92;
    m(2,0) = 11; m(2,1) = 69; m(2,2) =  5; m(2,3) = 86;
    m(3,0) =  8; m(3,1) =  9; m(3,2) = 98; m(3,3) = 23;

    Munkres<unsigned short> solver( m );
    std::cerr << solver();

    std::cerr << solver.cost( m ) << "\n";

    // TODO: Check...
    // 2,0
    // 1,1
    // 0,2
    // 3,3
    // Costs should be 140
  }

#endif
}
