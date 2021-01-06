#include "reader.h"

// explicit template member instantiation
template bool h5io::Reader::read( fields::Tile<1>& , ezh5::File& file);
template bool h5io::Reader::read( fields::Tile<2>& , ezh5::File& file);
template bool h5io::Reader::read( fields::Tile<3>& , ezh5::File& file);
template bool h5io::Reader::read( vlv::Tile<1>&    , ezh5::File& file);
template bool h5io::Reader::read( pic::Tile<2>&    , ezh5::File& file);
template bool h5io::Reader::read( pic::Tile<3>&    , ezh5::File& file);

