#ifndef ALEPH_UTILITIES_FILESYSTEM_HH__
#define ALEPH_UTILITIES_FILESYSTEM_HH__

#if defined(__unix__) || defined(__unix) || ( defined(__APPLE__) && defined(__MACH__) )
  #include <libgen.h>
  #include <unistd.h>

  #include <sys/stat.h>
#endif

#if defined(_BSD_SOURCE) || defined(_SVID_SOURCE) || defined(_DEFAULT_SOURCE)
  #define ALL_FILE_TYPE_QUERIES_AVAILABLE
#endif

#if defined(ALL_FILE_TYPE_QUERIES_AVAILABLE) || defined(_BSD_SOURCE) || defined(_DEFAULT_SOURCE) || ( defined(_XOPEN_SOURCE) && _XOPEN_SOURCE >= 500 )
  #define SOCKET_FILE_TYPE_QUERY_AVAILABLE
#endif

namespace aleph
{

namespace utilities
{

/**
  Describes potential file types. This is used *internally* but also by
  methods such as the `fileType()` function.
*/

enum class FileType
{
  BlockDevice,
  CharacterDevice,
  Directory,
  NamedPipe,
  RegularFile,
  Socket,
  SymbolicLink,
  Undefined
};

/** Returns the file type of a path */
FileType fileType( const std::string& path )
{
  FileType t = FileType::Undefined;

#if defined(__unix__) && ( ( defined(_POSIX_C_SOURCE) && _POSIX_C_SOURCE >= 200112L ) || defined(_DEFAULT_SOURCE) )
  struct stat info;
  int error = lstat( path.c_str(), &info );

  if( error == 0 )
  {
  #if defined(ALL_FILE_TYPE_QUERIES_AVAILABLE) || defined(_XOPEN_SOURCE)
    if( S_ISBLK( info.st_mode ) )
      t = FileType::BlockDevice;
    else if( S_ISCHR( info.st_mode ) )
      t = FileType::CharacterDevice;
    else if( S_ISDIR( info.st_mode ) )
      t = FileType::Directory;
    else if( S_ISFIFO( info.st_mode ) )
      t = FileType::NamedPipe;
    else if( S_ISREG( info.st_mode ) )
      t = FileType::RegularFile;
    else if( S_ISLNK( info.st_mode ) )
      t = FileType::SymbolicLink;
    #ifdef SOCKET_FILE_TYPE_QUERY_AVAILABLE
    else if( S_ISSOCK( info.st_mode ) )
      t = FileType::Socket;
    #endif
  #endif
  }
#endif

  return t;
}

/** Checks whether a given path is a directory */
bool isDirectory( const std::string& path )
{
  return fileType( path ) == FileType::Directory;
}

/** Checks whether a given path is a regular file */
bool isRegularFile( const std::string& path )
{
  return fileType( path ) == FileType::RegularFile;
}

/** Checks whether a given path is a socket */
bool isSocket( const std::string& path )
{
  return fileType( path ) == FileType::Socket;
}

/** Checks whether a path or a file exists */
bool exists( const std::string& path )
{
  return isDirectory( path ) || isRegularFile( path ) || isSocket( path );
}

/** Returns the basename, i.e the filename portion, of a path */
std::string basename( const std::string& path )
{
  if( path.empty() )
    return {};

  std::string result;

#if defined(_POSIX_VERSION) && _POSIX_VERSION >= 200112L

  char* buffer = new char[ path.size() + 1 ];
  std::copy( path.begin(), path.end(), buffer );

  buffer[ path.size() ] = '\0'; // `basename()` expects the string to be
                                // null-terminated

  char* name = ::basename( buffer );
  if( name )
    result = name;

  // `basename()` is guaranteed not to allocate _more_ memory. It might possibly
  // return a pointer to statically allocated memory (which we must not free) or
  // to parts of the input `path` pointer.
  //
  // Freeing said pointer is thus sufficient.
  delete[] buffer;

#else
  #error "No compatible implementation of 'basename' available"
#endif

  return result;
}

/**
  @returns Stem of the path. The stem is \i either the complete filename in
  case the path does not contain a dot \i or that part of the filename that
  precedes it. Hence, \c /foo/bar.txt will have a stem of \c bar.

  By definition, the special directories \c . and \c .. will remain their
  own stems. Normally, the stem will not contain a dot, though.
*/

std::string stem( const std::string& path )
{
  auto&& filename = basename( path );

  if( filename.find( '.' ) == std::string::npos || filename == "." || filename == ".." )
    return filename;

  auto pos = filename.find_last_of( '.' );
  return filename.substr( 0, pos );
}

/**
  @returns File extension of a path. This only works if the filename portion
  of the path contains a dot but is neither \c . nor \c .., i.e. the current
  directory or the parent directory. Note that the returned value will start
  with a dot in order to permit differentiating between filenames without an
  extension and filenames with an empty extension, i.e. a single dot. Though
  this behaviour may seem strange, it is compliant with \c Boost.Filesystem,
  for example.

  @see http://permalink.gmane.org/gmane.comp.lib.boost.devel/199744
*/

std::string extension( const std::string& path )
{
  auto&& filename = basename( path );

  if( filename.find( '.' ) == std::string::npos || filename == "." || filename == ".." )
    return {};

  auto pos = filename.find_last_of( '.' );
  return filename.substr( pos );
}

} // namespace utilities

} // namespace aleph

#endif
