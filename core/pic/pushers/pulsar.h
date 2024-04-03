#pragma once

#include "core/pic/pushers/pusher.h"

namespace pic {

/// Tailor-made pusher for pulsar setup. 
///
/// Based on a modified reduced guiding center approximation 
/// pusher (https://arxiv.org/abs/1701.05605) but with built-in
/// field interpolation and gravity
//
template<size_t D, size_t V>
class PulsarPusher :
  public Pusher<D,V>
{
  public:

  // configuration parameters
  double rad_star = 10.0;
  double rad_pcap = 1.0;
  double period_star = 0.0;

  // dipole parameters
  //double B0 = 1.0;         // Initial magnetic field strength B_0
  //double chi = 0.0;        // Obliquity angle between rotation axis and magnetic moment
  //double phase = 0.0;      // rotator phase
  //double delta_pc  = 1.0;  // polar smoothing

  double cenx, ceny, cenz; // center of the sphere

  double gravity_const = 1.0; // g_0 surface gravity constant controlling strength at the acceleration eq

  // pusher
  void push_container(
          pic::ParticleContainer<D>& container, 
          pic::Tile<D>& tile) override;
};

} // end of namespace pic
