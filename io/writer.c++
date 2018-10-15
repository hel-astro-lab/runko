#include "writer.h"

// explicit template member instantiation
template bool h5io::Writer::write(const fields::Tile<1>&, WriteMode::Standard);
template bool h5io::Writer::write(const fields::Tile<2>&, WriteMode::Standard);

template bool h5io::Writer::write(const fields::Tile<1>&, WriteMode::Analysis);
template bool h5io::Writer::write(const fields::Tile<2>&, WriteMode::Analysis);

template bool h5io::Writer::write(const vlv::Tile<1>&,    WriteMode::Standard);
