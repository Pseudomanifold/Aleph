import os
import re
import subprocess
import sys

from setuptools                   import setup
from setuptools                   import Extension
from setuptools.command.build_ext import build_ext
from distutils.version            import LooseVersion

class SourceExtension(Extension):
    """
    Describes a new extension that is built by source. This is required
    because the `Extension` base class assumes that an extension module
    consists of a *single* source file only.
    """
    def __init__(self, name):
        Extension.__init__(self,
                           name,
                           sources = [] ) # Do not supply any sources here
                                          # because everything is built by
                                          # CMake.

class Build(build_ext):
    """
    Class to enable an external build process to happen using CMake, making
    it possible to build the module directly without additional build steps
    required.
    """

    def run(self):
        try:
            subprocess.check_output(['cmake', '--version'])
        except subprocess.CalledProcessError:
            raise RuntimeError("CMake check failed")

        for extension in self.extensions:
            self.build_extension(extension)

    def build_extension(self, extension):
        extension_dir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(extension.name)))
        cmake_args    = ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=' + extension_dir,
                         '-DPYTHON_EXECUTABLE=' + sys.executable]

        cfg         = 'Debug' if self.debug else 'Release'
        #build_args  = ['--config', cfg]
        cmake_args += ['-DCMAKE_BUILD_TYPE=' + cfg]

        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)

        # TODO: fix path specification
        subprocess.check_call(['cmake', os.path.abspath('../../../../')] + cmake_args,
            cwd = self.build_temp,
            env = os.environ.copy()
        )

        subprocess.check_call(['cmake', '--build', '.', '--target', extension.name], #+ build_args,
            cwd = self.build_temp
        )

setup(
    name        = 'Aleph',
    version     = '${PACKAGE_VERSION}', # TODO: might want to use commit ID here
    packages    = [ 'aleph' ],
    package_dir = {
      '': '${CMAKE_CURRENT_BINARY_DIR}'
    },
    # FIXME: this does not actually build the package from source now
    # because the process clashes with `pipenv`
    ext_modules = [],
    cmdclass    = dict(build_ext=Build),
    zip_safe    = False
)
