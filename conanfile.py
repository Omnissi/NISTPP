from conans import ConanFile, CMake, tools
import version

class NistppConan(ConanFile):
    name = "NISTPP"
    license = "Unlicense"
    author = "Negodyaev Sergey (negodyaev.sergey@outlook.com)"
    description = "NIST test oc C++"
    settings = "os", "compiler", "build_type", "arch"
    options = {"shared": [True, False], "fPIC": [True, False]}
    default_options = {"shared": False, "fPIC": True}
    generators = ["cmake", "cmake_find_package", "cmake_paths"]
    exports = ["CMakeLists.txt", "version.py" "cmake/*"]

    def set_version(self):
        self.version = version.get_version()

    def build(self):
        cmake = CMake(self)
        cmake.configure()
        cmake.build()

    def package(self):
        self.copy("*.h", dst="include", keep_path=False)
        self.copy("*.hpp", dst="include", keep_path=False)
        self.copy("*.lib", dst="lib", keep_path=False)
        self.copy("*.dll", dst="bin", keep_path=False)
        self.copy("*.so", dst="lib", keep_path=False)
        self.copy("*.dylib", dst="lib", keep_path=False)
        self.copy("*.a", dst="lib", keep_path=False)