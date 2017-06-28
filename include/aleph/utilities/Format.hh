#ifndef ALEPH_UTILITIES_FORMAT_HH__
#define ALEPH_UTILITIES_FORMAT_HH__

#include <iomanip>
#include <ostream>
#include <string>

#include <cmath>

namespace aleph
{

namespace utilities
{

/**
  Formats a number according to an expected maximum number, padding
  unused digits with a fill character. For example, if the input is
  $n=5$ and the maximum is $N=999$, a fill character of `0` results
  in the an output of `005`.

  This utility function is useful for generating file names.
*/

template <class T> std::string format( T n, T N, char fill = '0' )
{
  int width = static_cast<int>(std::log10( N ) + 1);

  std::ostringstream stream;
  stream << std::setw( width ) << std::setfill( fill ) << n;

  return stream.str();
}

} // namespace utilities

} // namespace aleph

#endif
