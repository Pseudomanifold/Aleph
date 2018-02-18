#include <tests/Base.hh>

#include <aleph/utilities/Filesystem.hh>

#include <fstream>

void testUtilities()
{
  ALEPH_TEST_BEGIN( "Filesystem utilities" );

  std::string path1 = "/the/path/to/hell/is/paved/with/good/queries.txt";

  ALEPH_ASSERT_THROW( aleph::utilities::basename( path1 )  == std::string( "queries.txt" ) );
  ALEPH_ASSERT_THROW( aleph::utilities::stem( path1 )      == std::string( "queries" ) );
  ALEPH_ASSERT_THROW( aleph::utilities::extension( path1 ) == std::string( ".txt" ) );

  auto tmp = aleph::utilities::tempDirectory();

  ALEPH_ASSERT_THROW( tmp.empty() == false );

  ALEPH_TEST_END();
}

void testFileType()
{
  ALEPH_TEST_BEGIN( "File type detection" );

  auto tempDirectory = aleph::utilities::tempDirectory();
  auto test = tempDirectory + "/aleph_test.txt";

  {
    std::ofstream out( test );

    ALEPH_ASSERT_THROW( out );
  }

  ALEPH_ASSERT_THROW( aleph::utilities::exists( test ) );
  ALEPH_ASSERT_THROW( aleph::utilities::isRegularFile( test ) );

  ALEPH_TEST_END();
}

int main( int, char** )
{
  testUtilities();
  testFileType();
}
