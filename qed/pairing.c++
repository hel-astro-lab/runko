#include "pairing.h"

template class qed::Pairing<2>;
template class qed::Pairing<3>;


// explicit instantion of templated class method
//template void qed::Pairing::solve<2>(pic::Tile<2>& );
//template void qed::Pairing::solve<3>(pic::Tile<3>& );
//
//template void qed::Pairing::comp_pmax<2>(string, float_p, std::map<std::string, pic::ParticleContainer<2>*>& cons);
//template void qed::Pairing::comp_pmax<3>(string, float_p, std::map<std::string, pic::ParticleContainer<3>*>& cons);
//
//template void qed::Pairing::solve_mc<2>(pic::Tile<2>& );
//template void qed::Pairing::solve_mc<3>(pic::Tile<3>& );
//
//template void qed::Pairing::rescale<2>(pic::Tile<2>&, string&, double );
//template void qed::Pairing::rescale<3>(pic::Tile<3>&, string&, double );
