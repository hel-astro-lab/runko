#include "io/writers/writer.h"

// explicit template member instantiation
  
// emf
//template bool h5io::Writer::write(const emf::Tile<1>&, ezh5::File&); 
//template bool h5io::Writer::write(const emf::Tile<2>&, ezh5::File&); 
template bool h5io::Writer::write(const emf::Tile<3>&, ezh5::File&); 

// pic
//template bool h5io::Writer::write(const pic::Tile<1>&   , ezh5::File&); 
//template bool h5io::Writer::write(const pic::Tile<2>&   , ezh5::File&); 
template bool h5io::Writer::write(const pic::Tile<3>&   , ezh5::File&); 
