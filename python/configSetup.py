import numpy as np
from configparser import ConfigParser
import ast



section_names = 'io', 'simulation', 'grid', 'vmesh', 'particles', 'problem'


# Class for parsing and holding all configuration files
#
# see https://stackoverflow.com/questions/12620602/creating-a-class-with-all-the-elements-specified-in-a-file-using-configparser
# after this it is modified to automatically typecast into safe python datatypes
class Configuration(object):

    def __init__(self, *file_names):
        parser = ConfigParser()
        parser.optionxform = str  # make option names case sensitive

        found = parser.read(file_names)
        if not found:
            raise ValueError('No config file found!')
        for name in section_names:

            # here the magic happens
            elems = parser.items(name)
            for elem in elems:
                key = elem[0]
                val = ast.literal_eval(elem[1])

                self.__dict__.update( {key: val} )


        if not("dx" in self.__dict__):
            self.__dict__["dx"] = 1.0


        #compute dt from CFL (if given)
        if "cfl" in self.__dict__:
            self.__dict__["dt"] = self.__dict__["cfl"]*self.__dict__["dx"]
        #else:
        #    self.__dict__[key] = xx


        ##################################################
        # Default parameters; if not given

        if not("gamma_e" in self.__dict__):
            self.__dict__["gamma_e"] = 0.0
        if not("gamma_i" in self.__dict__):
            self.__dict__["gamma_i"] = 0.0


        #Compute initial values for configuration; 
        # assume that flow speed is beta if value is <=1, else it is gamma
        if np.abs(self.gamma_e) <= 1.0:
            self.ub_e = np.abs(self.gamma_e)
        else:
            self.ub_e = np.sqrt(self.gamma_e**2.0 - 1.0)
    
        if np.abs(self.gamma_i) <= 1.0:
            self.ub_i = np.abs(self.gamma_i)
        else:
            self.ub_i = np.sqrt(self.gamma_i**2.0 - 1.0)

        # imprint sign back
        self.ub_e *= np.sign(self.gamma_e)
        self.ub_i *= np.sign(self.gamma_i)

        #self.gamma_e = np.abs(self.gamma_e) 
        #self.gamma_i = np.abs(self.gamma_i) 

        # assume 1D/2D dimensionality if not otherwise specified
        if not("Ny" in self.__dict__):
            self.__dict__["Ny"] = 1

        if not("Nz" in self.__dict__):
            self.__dict__["Nz"] = 1

        #if velocity mesh cell size is only given, compute size
        if not("Nvx" in self.__dict__):
            self.__dict__["Nvx"] = int((self.vxmax - self.vxmin) // self.dvx)
        if not("Nvy" in self.__dict__):
            self.__dict__["Nvy"] = int((self.vymax - self.vymin) // self.dvy)
        if not("Nvz" in self.__dict__):
            self.__dict__["Nvz"] = int((self.vzmax - self.vzmin) // self.dvz)

        print("Nvx=", self.Nvx)



#load default 
#conf = Configuration('config.ini')



