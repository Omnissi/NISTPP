#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os

class GitVersion():
    def __init__(self):
        os.chdir(os.path.dirname(os.path.realpath(__file__)))

    def version(self):
        stream = os.popen("git describe --match [0-9]*.[0-9]*.[0-9]* --abbrev=0 --tags")
        return stream.read().strip()

    def build_number(self):
        stream = os.popen("git rev-list {}.. --count".format(self.version()))
        return stream.read().strip()

    def full_version(self):
        return "{}.{}".format(self.version(), self.build_number())
