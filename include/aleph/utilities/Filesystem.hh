#ifndef ALEPH_UTILITIES_FILESYSTEM_HH__
#define ALEPH_UTILITIES_FILESYSTEM_HH__

// If either one of these is defined, there is a good chance that POSIX
// concepts are available under the current architecture.
#if defined(__unix__) || defined(__unix) || ( defined(__APPLE__) && defined(__MACH__) )
  #define POSIX_SOURCES_AVAILABLE
#endif

#ifdef POSIX_SOURCES_AVAILABLE
  #include <libgen.h>
  #include <paths.h>
  #include <stdlib.h>
  #include <unistd.h>

  #include <sys/stat.h>
#endif

// In the best case, the `_POSIX_VERSION` variable is set on the system,
// but more exotic configurations may not have this. The essence of this
// check is to ensure that common features of `lstat()` are available.
#if ( defined(_POSIX_VERSION) && _POSIX_VERSION >= 200112L ) || ( defined(POSIX_SOURCES_AVAILABLE) && ( ( defined(_POSIX_C_SOURCE) && _POSIX_C_SOURCE >= 200112L ) || defined(_DEFAULT_SOURCE) ) )
  #define POSIX_VERSION_COMPATIBLE
#endif

#if defined(_BSD_SOURCE) || defined(_SVID_SOURCE) || defined(_DEFAULT_SOURCE)
  #define ALL_FILE_TYPE_QUERIES_AVAILABLE
#endif

// Brute-force check for the availability of some macros for `lstat()`,
// in case the check from above fails.
#if defined(S_ISBLK) && defined(S_ISCHR) && defined(S_ISDIR) && defined(S_ISFIFO) && defined(S_ISREG) && defined(S_ISLNK) && defined(S_ISSOCK)
  #define ALL_FILE_TYPE_QUERIES_AVAILABLE
#endif

#if defined(ALL_FILE_TYPE_QUERIES_AVAILABLE) || defined(_BSD_SOURCE) || defined(_DEFAULT_SOURCE) || ( defined(_XOPEN_SOURCE) && _XOPEN_SOURCE >= 500 )
  #define SOCKET_FILE_TYPE_QUERY_AVAILABLE
#endif

#include <fstream>

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

#ifdef POSIX_VERSION_COMPATIBLE
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
  return    isDirectory( path )
         || isRegularFile( path )
         || isSocket( path )
         || ( std::ifstream( path ) ); // fall-back in case no other queries are available;
                                       // we just try to open the file instead
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

/**
  @returns The temporary directory of the given system. This function
  attempts to solve this in a platform-independent manner. It is very
  likely though that this only works for POSIX-based systems.
*/

std::string tempDirectory()
{
#ifdef POSIX_SOURCES_AVAILABLE
  auto tmpdir  = getenv("TMPDIR");
  auto tmp     = getenv("TMP");
  auto temp    = getenv("TEMP");
  auto tempdir = getenv("TEMPDIR");

  // These may not always be available, so we have to perform an
  // additional check below.
  const char* ptmpdir = nullptr;
  const char* pathtmp = nullptr;

  #ifdef P_tmpdir
    ptmpdir = P_tmpdir;
  #endif

  #ifdef _PATH_TMP
    pathtmp = _PATH_TMP;
  #endif

  // This does not check to what extent the information present in two
  // of these variables is consistent with the rest. The order is more
  // or less random here; I was unable to find an authoritative guide.

  return tmpdir ? tmpdir
                : tmp ? tmp
                      : temp ? temp
                             : tempdir ? tempdir
                                       : ptmpdir ? ptmpdir
                                                 : pathtmp ? pathtmp
                                                           : ""; // Return an empty directory rather
                                                                 // than an incorrect one
#endif

  // Return an empty directory rather than guessing an incorrect one;
  // we *could* potentially default to '/tmp' here, but I do not like
  // this idea of an untested fallback.
  return {};
}

} // namespace utilities

} // namespace aleph

#endif
