#ifndef COMMON_H
#define COMMON_H

#include "definitions.hpp"


// physical stuff
//-------------------------------------------------- 
const double c = 10.0; ///< Speed of light
const double q = 1.0; ///< charge of ions
const double e = 1.0;  ///< elementary charge
const double pi = 3.14159265359;  ///<  mathematical pi


const double me = 1.0; ///<  mass of electron
const double mp = 16.0; ///<  mass of proton


const std::vector<CellID>& getLocalCells();
void recalculateLocalCellsCache();



#endif
