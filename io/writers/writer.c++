#include "writer.h"

// explicit template member instantiation
  
// fields
template bool h5io::Writer::write(const fields::Tile<1>&, ezh5::File&); 
template bool h5io::Writer::write(const fields::Tile<2>&, ezh5::File&); 
template bool h5io::Writer::write(const fields::Tile<3>&, ezh5::File&); 

// vlasov
template bool h5io::Writer::write(const vlv::Tile<1>&   , ezh5::File&); 

// pic
template bool h5io::Writer::write(const pic::Tile<2>&   , ezh5::File&); 
template bool h5io::Writer::write(const pic::Tile<3>&   , ezh5::File&); 
