#!/usr/bin/env python
# -*- coding: utf-8 -*-
from conans import ConanFile
from conans.tools import load
import re

def get_version():
    try:
        content = load("CMakeLists.txt")
        version = re.search(r"CXX VERSION (.*)\)", content).group(1)
        return version.strip()
    except Exception as e:
        return None
