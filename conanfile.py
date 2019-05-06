from conans import ConanFile, CMake, tools


class AlephConan(ConanFile):
    name = "Aleph"
    version = "0.1.0"
    license = "MIT"
    author = "Bastian Rieck <bastian@rieck.ru>"
    url = "https://github.com/Pseudomanifold/Aleph"
    description = "A library for exploring persistent homology"
    topics = ("persistent-homology", "topological-data-analysis", "tda")
    settings = "os", "compiler", "build_type", "arch"
    options = {
        "enable_rapidjson": [True, False],
        "enable_pybind11": [True, False],
        "enable_eigen3": [True, False],
        "enable_tinyxml2": [True, False]
    }

    default_options = {
        "enable_rapidjson=False",
        "enable_pybind11=False",
        "enable_eigen3=False",
        "enable_tinyxml2=False"
        }

    generators = "cmake"

    def requirements(self):
        self.requires("boost/[>1.55]@conan/stable")

        if self.options.enable_rapidjson:
            self.requires("rapidjson/[>1.0]@bincrafters/stable")

        # FIXME: version is not optimal
        if self.options.enable_pybind:
            self.requires("pybind11/2.2.4@conan/stable")

        # FIXME: version is not optimal
        if self.options.enable_eigen:
            self.requires("eigen/3.3.5@conan/stable")

        # FIXME: version is not optimal
        if self.options.enable_tinyxml2:
            self.requires("tinyxml2/7.0.1@nicolastagliani/stable")

    def source(self):
        self.run("git clone https://github.com/Pseudomanifold/Aleph.git")
        self.run("cd Aleph")

    def build(self):
        cmake = CMake(self)
        cmake.configure(source_folder="Aleph")
        cmake.build()

    def package(self):
        self.copy("*.hh", dst="include", src="Aleph")
        self.copy("*.so", dst="lib", keep_path=False)
        self.copy("*.dylib", dst="lib", keep_path=False)
        self.copy("*.a", dst="lib", keep_path=False)

    def package_info(self):
        self.cpp_info.libs = ["hello"]

