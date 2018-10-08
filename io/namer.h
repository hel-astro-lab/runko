#pragma once
#include <tuple>
#include <string>


namespace h5io {

using std::string;
using std::to_string;


/// General class to do all the text and file name mangling 
class Namer {

  private:
    const string extension = ".h5";


  public:
    string name;

    Namer(string& prefix) 
    {
      name = prefix + extension;
    }

    Namer(string& prefix, int lap) 
    {
      name = prefix + "_" + to_string(lap) + extension;
    }
      
    //standard numbering scheme
    // TODO generalize to variable argument
    string numbering(size_t i, size_t j, size_t k) const
    {
      return to_string(i) + "_" + to_string(j) + "_" + to_string(k);
    }

    string numbering(std::tuple<size_t, size_t, size_t>& ind) const
    {
      return numbering(
          std::get<0>(ind),
          std::get<1>(ind),
          std::get<2>(ind)
          );
    }

};


} // end of namespace io
