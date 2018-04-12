from distutils.core import setup

setup(
  name        = 'Aleph',
  version     = '${PACKAGE_VERSION}', # TODO: might want to use commit ID here
  packages    = [ 'aleph' ]
  package_dir = {
    '': '${CMAKE_CURRENT_BINARY_DIR}'
  },
  package_data = {
    '': ['aleph.so']
  }
)
