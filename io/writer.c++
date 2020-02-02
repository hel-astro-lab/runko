#include "writer.h"

// explicit template member instantiation
template bool h5io::Writer::write(const fields::Tile<1>& ); //, Write_mode::Standard);
template bool h5io::Writer::write(const fields::Tile<2>& ); //, Write_mode::Standard);

//template bool h5io::Writer::write2(const fields::Tile<1>& ); //, Write_mode::Analysis);
//template bool h5io::Writer::write2(const fields::Tile<2>& ); //, Write_mode::Analysis);

template bool h5io::Writer::write(const vlv::Tile<1>&    ); //,    Write_mode::Standard);

template bool h5io::Writer::write(const pic::Tile<2>&    ); //,    Write_mode::Standard);
