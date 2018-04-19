#include <aleph/persistenceDiagrams/PersistenceDiagram.hh>
#include <aleph/persistenceDiagrams/io/Raw.hh>

#include <string>
#include <vector>

#include <getopt.h>

using DataType           = double;
using PersistenceDiagram = aleph::PersistenceDiagram<DataType>;

int main( int argc, char** argv)
{
  bool normalize = false;

  {
    static option commandLineOptions[] =
    {
      { "normalize", no_argument, nullptr, 'n' },
      { nullptr    , 0          , nullptr,  0  }
    };

    int option = 0;

    while( ( option = getopt_long( argc, argv, "n", commandLineOptions, nullptr ) ) != -1 )
    {
      switch( option )
      {
      case 'n':
        normalize = true;
        break;
      default:
        break;
      }
    }
  }
}
