#include "writer.h"

// explicit template member instantiation
  
// fields
template bool h5io::Writer::write(const fields::Tile<1>& ); 
template bool h5io::Writer::write(const fields::Tile<2>& ); 
template bool h5io::Writer::write(const fields::Tile<3>& ); 

// vlasov
template bool h5io::Writer::write(const vlv::Tile<1>&    ); 

// pic
template bool h5io::Writer::write(const pic::Tile<2>&    ); 
template bool h5io::Writer::write(const pic::Tile<3>&    ); 
