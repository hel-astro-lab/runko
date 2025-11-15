#include "core/emf/yee_lattice.h"
#include "tyvi/mdgrid.h"

#include <tuple>

namespace emf {

void
  YeeLattice::push_b_FDTD2(const value_type dt)
{
  tyvi::mdgrid_work w {};
  push_b_FDTD2(w, dt);
  w.wait();
}

void
  YeeLattice::push_b_FDTD2(const tyvi::mdgrid_work& w, const value_type dt)
{
  /* FIXME: figure out if dt and corr from emf::FDTD2 are needed. */
  const auto h   = halo_size_;
  const auto hp1 = h + 1uz;

  const auto i_nonhalo = std::tuple { h, h + extents_wout_halo_[0] };
  const auto j_nonhalo = std::tuple { h, h + extents_wout_halo_[1] };
  const auto k_nonhalo = std::tuple { h, h + extents_wout_halo_[2] };

  const auto i_nonhalo_p1 = std::tuple { hp1, hp1 + extents_wout_halo_[0] };
  const auto j_nonhalo_p1 = std::tuple { hp1, hp1 + extents_wout_halo_[1] };
  const auto k_nonhalo_p1 = std::tuple { hp1, hp1 + extents_wout_halo_[2] };

  const auto Emds      = nonhalo_submds(E_.mds());
  const auto Bmds      = nonhalo_submds(B_.mds());
  const auto Emds_halo = E_.mds();

  const auto Emds_ip1 = std::submdspan(Emds_halo, i_nonhalo_p1, j_nonhalo, k_nonhalo);
  const auto Emds_jp1 = std::submdspan(Emds_halo, i_nonhalo, j_nonhalo_p1, k_nonhalo);
  const auto Emds_kp1 = std::submdspan(Emds_halo, i_nonhalo, j_nonhalo, k_nonhalo_p1);


  w.for_each_index(Bmds, [=](const auto idx) {
    const auto DkEy = Emds_kp1[idx][1] - Emds[idx][1];
    const auto DjEz = Emds_jp1[idx][2] - Emds[idx][2];
    Bmds[idx][0]    = Bmds[idx][0] + dt * (DkEy - DjEz);
  });
  w.for_each_index(Bmds, [=](const auto idx) {
    const auto DiEz = Emds_ip1[idx][2] - Emds[idx][2];
    const auto DkEx = Emds_kp1[idx][0] - Emds[idx][0];
    Bmds[idx][1]    = Bmds[idx][1] + dt * (DiEz - DkEx);
  });
  w.for_each_index(Bmds, [=](const auto idx) {
    const auto DjEx = Emds_jp1[idx][0] - Emds[idx][0];
    const auto DiEy = Emds_ip1[idx][1] - Emds[idx][1];
    Bmds[idx][2]    = Bmds[idx][2] + dt * (DjEx - DiEy);
  });
}

void
  YeeLattice::push_e_FDTD2(const value_type dt)
{
  tyvi::mdgrid_work w {};
  push_e_FDTD2(w, dt);
  w.wait();
}

void
  YeeLattice::push_e_FDTD2(const tyvi::mdgrid_work& w, const value_type dt)
{
  /* FIXME: figure out if dt and corr from emf::FDTD2 are needed. */

  const auto h   = halo_size_;
  const auto hm1 = static_cast<std::size_t>(halo_size_ - 1uz);

  const auto i_nonhalo = std::tuple { halo_size_, halo_size_ + extents_wout_halo_[0] };
  const auto j_nonhalo = std::tuple { halo_size_, halo_size_ + extents_wout_halo_[1] };
  const auto k_nonhalo = std::tuple { halo_size_, halo_size_ + extents_wout_halo_[2] };

  const auto i_nonhalo_m1 = std::tuple { hm1, hm1 + extents_wout_halo_[0] };
  const auto j_nonhalo_m1 = std::tuple { hm1, hm1 + extents_wout_halo_[1] };
  const auto k_nonhalo_m1 = std::tuple { hm1, hm1 + extents_wout_halo_[2] };

  const auto Emds      = nonhalo_submds(E_.mds());
  const auto Bmds      = nonhalo_submds(B_.mds());
  const auto Bmds_halo = B_.mds();

  const auto Bmds_im1 = std::submdspan(Bmds_halo, i_nonhalo_m1, j_nonhalo, k_nonhalo);
  const auto Bmds_jm1 = std::submdspan(Bmds_halo, i_nonhalo, j_nonhalo_m1, k_nonhalo);
  const auto Bmds_km1 = std::submdspan(Bmds_halo, i_nonhalo, j_nonhalo, k_nonhalo_m1);

  w.for_each_index(Emds, [=](const auto idx) {
    const auto DkBy = Bmds_km1[idx][1] - Bmds[idx][1];
    const auto DjBz = Bmds_jm1[idx][2] - Bmds[idx][2];
    Emds[idx][0]    = Emds[idx][0] + dt * (DkBy - DjBz);
  });
  w.for_each_index(Emds, [=](const auto idx) {
    const auto DiBz = Bmds_im1[idx][2] - Bmds[idx][2];
    const auto DkBx = Bmds_km1[idx][0] - Bmds[idx][0];
    Emds[idx][1]    = Emds[idx][1] + dt * (DiBz - DkBx);
  });
  w.for_each_index(Emds, [=](const auto idx) {
    const auto DjBx = Bmds_jm1[idx][0] - Bmds[idx][0];
    const auto DiBy = Bmds_im1[idx][1] - Bmds[idx][1];
    Emds[idx][2]    = Emds[idx][2] + dt * (DjBx - DiBy);
  });
}

}  // namespace emf
