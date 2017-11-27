#include <aleph/math/Bootstrap.hh>

#include <aleph/math/PiecewiseLinearFunction.hh>

#include <string>

#include <getopt.h>

int main( int argc, char** argv )
{
  double alpha = 0.95;

  {
    static option commandLineOptions[] =
    {
      { "alpha", required_argument, nullptr, 'a' },
      { nullptr, 0                , nullptr,  0  }
    };

    int option = 0;
    while( ( option = getopt_long( argc, argv, "a:", commandLineOptions, nullptr ) ) != - 1 )
    {
      switch( option )
      {
      case 'a':
        alpha = std::stod( optarg );
        break;
      }
    }
  }
}
