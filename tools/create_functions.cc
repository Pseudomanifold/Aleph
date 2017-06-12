#include <algorithm>
#include <iostream>
#include <numeric>
#include <random>
#include <vector>

int main(int, char**)
{
  unsigned min =  0;
  unsigned max = 20;

  std::random_device rd;
  std::mt19937 rng( rd() );

  for( unsigned k = 0; k < 500; k++ )
  {
    std::vector<unsigned> valuesIn( max );
    std::iota( valuesIn.begin(), valuesIn.end(), min+1 );
    std::shuffle( valuesIn.begin(), valuesIn.end(), rng );

    std::vector<unsigned> valuesOut = { 0 };

    bool previousIsMinimum = true;
    for( unsigned i = 0; i < 20; i++ )
    {
      auto previous = valuesOut.back();
      auto current  = previous;

      if( previousIsMinimum )
      {
        for( auto&& v : valuesIn )
        {
          if( v > previous )
          {
            current = v;
            break;
          }
        }
      }
      else
      {
        for( auto&& v : valuesIn )
        {
          if( v < previous )
          {
            current = v;
            break;
          }
        }
      }

      valuesIn.erase( std::remove( valuesIn.begin(), valuesIn.end(), current ),
                      valuesIn.end() );

      if( current != previous )
      {
        valuesOut.push_back( current );
        previousIsMinimum = !previousIsMinimum;
      }
    }

    valuesOut.push_back(0);

    std::cerr << "* Created " << valuesOut.size() << " function values\n";

    for( auto&& value : valuesOut )
      std::cout << value << " ";
    std::cout << "\n";
  }
}
