#ifndef ALEPH_TESTS_BASE_HH__
#define ALEPH_TESTS_BASE_HH__

#include <iostream>
#include <stdexcept>
#include <string>

namespace aleph
{

#define ALEPH_ASSERT_THROW( condition )                             \
{                                                                   \
  if( !( condition ) )                                              \
  {                                                                 \
    throw std::runtime_error(   std::string( __FILE__ )             \
                              + std::string( ":" )                  \
                              + std::to_string( __LINE__ )          \
                              + std::string( " in " )               \
                              + std::string( __PRETTY_FUNCTION__ )  \
    );                                                              \
  }                                                                 \
}

#define ALEPH_TEST_BEGIN( name )\
{\
  std::cerr << "-- Running test \"" << name << "\"...";\
}

#define ALEPH_TEST_END \
{\
  std::cerr << "finished\n";\
}


}

#endif
