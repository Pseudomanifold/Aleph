#include <iostream>
#include <string>

#include "distances/detail/Matrix.hh"
#include "distances/detail/Munkres.hh"

using namespace aleph::distances::detail;

int main()
{
  {
    Matrix<unsigned short> m( 3 );

    m(0,0) = 1; m(0,1) = 2; m(0,2) = 3;
    m(1,0) = 2; m(1,1) = 4; m(1,2) = 6;
    m(2,0) = 3; m(2,1) = 6; m(2,2) = 9;

    auto separator = std::string( 80, '-' );

    std::cerr << separator << "\n"
              << "Original costs\n"
              << separator << "\n"
              << m
              << separator << "\n";

    Munkres<unsigned short> solver( m );
    solver();

    std::cerr << separator << "\n"
              << "Modified costs\n"
              << separator << "\n"
              << m
              << separator << "\n";
  }

  {
    Matrix<unsigned short> m( 4 );

    m(0,0) = 82; m(0,1) = 83; m(0,2) = 69; m(0,3) = 92;
    m(1,0) = 77; m(1,1) = 37; m(1,2) = 49; m(1,3) = 92;
    m(2,0) = 11; m(2,1) = 69; m(2,2) =  5; m(2,3) = 86;
    m(3,0) =  8; m(3,1) =  9; m(3,2) = 98; m(3,3) = 23;

    Munkres<unsigned short> solver( m );
    std::cerr << solver();

    // TODO: Check...
    // 2,0
    // 1,1
    // 0,2
    // 3,3
    // Costs should be 140
  }
}
