from conans import ConanFile, CMake, tools
import version
import os

class NistppConan(ConanFile):
    name = "NISTPP"
    license = "Unlicense"
    author = "Negodyaev Sergey (negodyaev.sergey@outlook.com)"
    description = "NIST test oc C++"
    settings = "os", "compiler", "build_type", "arch"
    options = {"shared": [True, False], "fPIC": [True, False], "enable_tests": [True, False]}
    default_options = {"shared": False, "fPIC": True, "enable_tests": True}
    generators = ["cmake", "cmake_find_package", "cmake_paths"]
    exports = ["CMakeLists.txt", "version.py" "cmake/*"]

    def requirements(self):
        self.requires.add("boost/1.73.0")
        self.requires.add("eigen/3.4.0")
        if self.options.enable_tests:
            self.requires.add("gtest/cci.20210126")

    def set_version(self):
        self.version = version.get_version(os.path.dirname(os.path.abspath(__file__)))

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
