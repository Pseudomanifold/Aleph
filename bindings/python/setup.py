from distutils.core import setup

import sys
if sys.version_info < (3,0):
  sys.exit('Sorry, Python < 3.0 is not supported')

setup(
  name        = 'Aleph',
  version     = '${PACKAGE_VERSION}', # TODO: might want to use commit ID here
  packages    = [ 'aleph' ],
  package_dir = {
    '': '${CMAKE_CURRENT_BINARY_DIR}'
  },
  package_data = {
    '': ['aleph.so']
  },
  python_requires = '>3'
)
