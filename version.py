#!/usr/bin/env python
# -*- coding: utf-8 -*-
from conans import ConanFile
from conans.tools import load
import re

def get_version(source_dir: str):
    try:
        print(source_dir)
        content = load(source_dir + str("/CMakeLists.txt"))
        version = re.search(r"NISTPP VERSION (.*)\)", content).group(1)
        return version.strip()
    except Exception as e:
        print(e)
        return None
