# Module: check_compiler_features
# Author: Bastian Rieck <bastian.rieck@bsse.ethz.ch>
#
# Checks for some compiler features and sets variables accordingly. This
# is useful when building with:
#
# - different compilers (`gcc`, `clang`)
# - older versions of some compilers in an OS (looking at you, Ubuntu)
#
# Parts of these checks have been inspired by the great deal.II library
# that you can check out on GitHub (https://github.com/dealii/dealii).

#
# Checks whether we have access to the _Pragma directive and push/pop
# operations. If this is the case, we are able to disable some
# diagnostics for some external headers which decreases the amount of
# warnings.
#
CHECK_CXX_SOURCE_COMPILES(
  "
  _Pragma(\"GCC diagnostic push\")
  _Pragma(\"GCC diagnostic ignored \\\\\\\"-Wextra\\\\\\\"\")
  int main( int, char** )
  {
    return 0;
  }
  _Pragma(\"GCC diagnostic pop\")
  "
  ALEPH_COMPILER_HAS_DIAGNOSTIC_PRAGMA )

#
# Checks whether the regex token iterator is available. This is
# a feature of C++11 but not all compilers have implemented it.
# Some older versions of `gcc` only support parts of the C++11
# standard.
#

CHECK_CXX_SOURCE_COMPILES(
  "
  \#include <regex>

  int main( int, char** )
  {
    std::sregex_token_iterator t;
    return 0;
  }
  "
  ALEPH_COMPILER_HAS_REGEX_TOKEN_ITERATOR )

#
# Checks whether the compiler permits us to deprecate a function. If so,
# we can set up a macro that does this job for us. If not, this macro
# will be left empty.
#
CHECK_CXX_SOURCE_COMPILES(
  "
  int f() __attribute__((deprecated));
  int f()
  {
    return 42;
  }

  int main( int, char** )
  {
    return f();
  }
  "
  ALEPH_COMPILER_HAS_ATTRIBUTE_DEPRECATED )

IF( ALEPH_COMPILER_HAS_ATTRIBUTE_DEPRECATED )
  SET( ALEPH_DEPRECATED "__attribute__((deprecated))" )
ELSE()
  # This needs to contain at least a single space. Else, CMake will not let us
  # define this as a macro.
  SET( ALEPH_DEPRECATED " " )
ENDIF()
