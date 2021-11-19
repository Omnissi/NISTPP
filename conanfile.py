from conans import ConanFile, CMake, tools
from version import GitVersion
import os

class NistppConan(ConanFile):
    name = "nistpp"
    license = "Unlicense"
    author = "Negodyaev Sergey (negodyaev.sergey@outlook.com)"
    description = "NIST test on C++"
    url = "https://git.omnissi-factory.ru/Omnissi/NISTPP"
    settings = "os", "compiler", "build_type", "arch", "os_build"
    options = {"shared": [True, False], "fPIC": [True, False], "enable_tests": [True, False]}
    default_options = {"shared": False, "fPIC": True, "enable_tests": True}
    generators = ["cmake", "cmake_find_package", "cmake_paths"]

    exports_sources = ["nistpp/*", "Sprout/*", "unit_tests/*", "cmake/*", "CMakeLists.txt"]
    exports = ["conanfile.py", "version.py", "version.h.in"]
    _cmake = None

    def requirements(self):
        self.requires.add("boost/1.73.0")
        self.requires.add("kissfft/131.1.0")
        if self.options.enable_tests:
            self.requires.add("gtest/cci.20210126")

    def configure(self):
        self.options["kissfft"].openmp = True
        self.options["kissfft"].datatype = "double"

    def _configure_cmake(self):
        if self._cmake:
            return self._cmake
        cmake = CMake(self)
        cmake.definitions["ENABLE_TEST"] = self.options.enable_tests
        cmake.definitions["OS_BUILD"] = self.settings.os_build
        cmake.configure()
        self._cmake = cmake
        return self._cmake

    def set_version(self):
        self.version = GitVersion().full_version()

    def build(self):
        cmake = self._configure_cmake()
        cmake.build()

    def package(self):
        cmake = self._configure_cmake()
        cmake.install()

    def package_info(self):
        self.cpp_info.build_modules = [ "cmake/nistpp.cmake" ]
        self.cpp_info.libs = ["nistpp"]