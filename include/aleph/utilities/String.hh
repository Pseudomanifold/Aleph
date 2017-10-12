#ifndef ALEPH_UTILITIES_STRING_HH__
#define ALEPH_UTILITIES_STRING_HH__

#include <algorithm>
#include <limits>
#include <regex>
#include <sstream>
#include <string>
#include <vector>

#include <cctype>

namespace aleph
{

namespace utilities
{

template <class T> T ltrim( T sequence )
{
  sequence.erase( sequence.begin(),
                  std::find_if_not( sequence.begin(), sequence.end(),
                                    [] ( const typename T::value_type& c )
                                    {
                                      return std::isspace( c );
                                    } ) );

  return sequence;
}

template <class T> T rtrim( T sequence )
{
  sequence.erase( std::find_if_not( sequence.rbegin(), sequence.rend(),
                                    [] ( const typename T::value_type& c )
                                    {
                                      return std::isspace( c );
                                    } ).base(),
                  sequence.end() );

  return sequence;
}

template <class T> T trim( T sequence )
{
  return ltrim( rtrim( sequence ) );
}

template <class T> std::vector<T> split( const T& sequence,
                                         const T& regex = "[[:space:]]+" )
{
  std::regex re( regex );
  std::sregex_token_iterator begin( sequence.begin(), sequence.end(), re, -1 );
  std::sregex_token_iterator end;

  return { begin, end };
}

/**
  Attempts to convert a sequence type `S` to a non-sequence type `T` by
  using `std::stringstream`. This makes converting strings to different
  types such as numbers easier.

  @tparam S Sequence type (e.g. `std::string`)
  @tparam T Non-sequence type (e.g. `ìnt`)

  @param sequence Sequence to convert
  @param success  Flag indicating the success of the conversion

  @returns Result of the conversion. Errors do *not* result in an error
  being thrown. Use the \p success parameter to check for errors.
*/

template <class T, class S> T convert( const S& sequence, bool& success )
{
  T result = T();
  success  = false;

  std::istringstream converter( sequence );
  converter >> result;

  // Try some special handling for some special tokens. Other errors are
  // silently ignored. I am not sure whether this is the right behaviour
  // but I see no pressing reason to change it now.
  if( converter.fail() )
  {
    std::string string = sequence;

    std::transform( string.begin(), string.end(),
                    string.begin(), ::tolower );

    if( string == "+inf" || string == "inf" || string == "+infinity" || string == "infinity" )
      result = std::numeric_limits<T>::infinity();
    else if ( string == "-inf" || string == "-infinity" )
      result = -std::numeric_limits<T>::infinity();
    else if( string == "nan" )
      result = std::numeric_limits<T>::quiet_NaN();

    success = result != T();
  }
  else
    success = true;

  return result;

}

/** @overload convert() */
template <class T, class S> T convert( const S& sequence )
{
  bool success = false;
  return convert<T>( sequence, success );
}

} // namespace utilities

} // namespace aleph

#endif
