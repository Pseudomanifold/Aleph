#include <tests/Base.hh>

#include <aleph/utilities/Filesystem.hh>

int main( int, char** )
{
  ALEPH_TEST_BEGIN( "Filesystem utilities" );

  std::string path1 = "/the/path/to/hell/is/paved/with/good/queries.txt";

  ALEPH_ASSERT_THROW( aleph::utilities::basename( path1 )  == std::string( "queries.txt" ) );
  ALEPH_ASSERT_THROW( aleph::utilities::stem( path1 )      == std::string( "queries" ) );
  ALEPH_ASSERT_THROW( aleph::utilities::extension( path1 ) == std::string( ".txt" ) );

  ALEPH_TEST_END();
}
