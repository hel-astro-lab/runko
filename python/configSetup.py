from ConfigParser import SafeConfigParser
import ast


section_names = 'io', 'simulation', 'grid', 'vmesh', 'particles'


# Class for parsing and holding all configuration files
#
# see https://stackoverflow.com/questions/12620602/creating-a-class-with-all-the-elements-specified-in-a-file-using-configparser
# after this it is modified to automatically typecast into safe python datatypes
class Configuration(object):

    def __init__(self, *file_names):
        parser = SafeConfigParser()
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
        
        #compute dt from CFL (if given)
        if "cfl" in self.__dict__:
            self.__dict__["dt"] = self.__dict__["cfl"]*self.__dict__["dx"]
        #else:
        #    self.__dict__[key] = xx




#load default 
#conf = Configuration('config.ini')



