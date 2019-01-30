from conans import ConanFile, CMake, tools

class AlephConan(ConanFile):
    name = "Aleph"
    version = "0.0.1"
    license = "MIT"
    url = "https://github.com/Submanifold/Aleph"
    description = "A library for exploring persistent homology https://submanifold.github.io/Aleph/"
    settings = "os", "compiler", "build_type", "arch"
    options = {"enable_rapidjson" : [True, False], "enable_pybind" : [True, False], "enable_eigen" : [True, False], "enable_tinyxml2" : [True, False]}
    default_options = "enable_rapidjson=False", "enable_pybind=False", "enable_pybind=False", "enable_eigen=False", "enable_tinyxml2=False" 
    generators = "cmake"
    exports_sources = ["Findtinyxml2.cmake"]

    def requirements(self):
        self.requires("boost/1.65.1@conan/stable")

        if self.options.enable_rapidjson:
            self.requires("rapidjson/1.1.0@bincrafters/stable")

        if self.options.enable_pybind:
            self.requires("pybind11/2.2.4@conan/stable")

        if self.options.enable_eigen:
            self.requires("eigen/3.3.5@conan/stable")

        if self.options.enable_tinyxml2:
            self.requires("tinyxml2/7.0.1@nicolastagliani/stable")

    def package_info(self):
        self.copy("Findtinyxml2.cmake", ".", ".")
