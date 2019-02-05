#pragma once

#include <vector>
#include <string>

#include "namer.h"
#include "../corgi/corgi.h"


namespace h5io { 


template<size_t D>
class QuickWriter {

  private:

    /// general file extension to be appended to file names
    const string extension = ".h5";
    
    /// Object to handle file names and extensions
    std::string fname;

  public:

    /// constructor that creates a name and opens the file handle
    QuickWriter(const std::string& prefix, int lap) 
    {
      fname = prefix + "-" + to_string(lap) + extension;
    }

    bool write(corgi::Node<D>& grid);
};


} // end of namespace h5io

