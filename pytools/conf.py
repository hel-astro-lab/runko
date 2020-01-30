# -*- coding: utf-8 -*-

from configparser import ConfigParser
import ast

section_names = "io", "simulation", "grid", "vmesh", "particles", "problem"


# Class for parsing and holding all configuration files
#
# see https://stackoverflow.com/questions/12620602/creating-a-class-with-all-the-elements-specified-in-a-file-using-configparser
# after this it is modified to automatically typecast into safe python datatypes
class Configuration(object):
    def __init__(self, *file_names, **kwargs):

        # check that we do not get an empty file name
        for file_name in file_names:
            if file_name == None:
                raise Exception("No configuration file (--conf) provided")

        # built parser and read .ini
        parser = ConfigParser()
        parser.optionxform = str  # make option names case sensitive

        found = parser.read(file_names)

        if not found:
            raise ValueError("No config file found!")
        for name in section_names:

            # here the magic happens and we turn parameters into conf members
            elems = parser.items(name)
            for elem in elems:
                key = elem[0]
                val = ast.literal_eval(elem[1])
                self.__dict__.update({key: val})

        # automatically calculate grid size
        self.xmin = 0.0
        self.ymin = 0.0
        self.zmin = 0.0
        self.xmax = self.Nx * self.NxMesh
        self.ymax = self.Ny * self.NyMesh
        self.zmax = self.Nz * self.NzMesh

        # dimensionality switches
        self.oneD = False
        self.twoD = False
        self.threeD = False

        if self.Nz > 1:
            self.threeD = True
        elif self.Ny > 1:
            self.twoD = True
        elif self.Nx > 1:
            self.oneD = True
