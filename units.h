#pragma once

#include <cmath>

/*! See https://github.com/awsteiner/o2scl/src/base/constants.h
 *
 * CODATA 2010 values are in Mohr12. IAU 2009 values are from Luzum11. 
 * Solar mass from http://asa.usno.navy.mil/SecK/2013/Astronomical_Constants_2013.pdf
 *
 */



struct cgs {
  public:
    /// \f$ \pi \f$ 
    const double pi=std::acos(-1.0);
    /// cm / s
    const double speed_of_light=2.99792458e10;
    /// Newtonian constant of gravitation in cm^3 / g s^2 (CODATA 2010 value)
    const double gravitational_constant=6.67384e-8;
    /// Planck constant in g cm^2 / s (CODATA 2010 value)
    const double plancks_constant_h=6.62606957e-27;
    /// Astronomical unit in cm (IAU 2009 value; now exact)
    const double astronomical_unit=1.49597870700e13;
    /// cm
    const double light_year=9.46053620707e17;
    /// cm
    const double parsec=3.08567758135e18;
    /// cm / s^2
    const double grav_accel=9.80665e2;
    /// Electron volt in g cm^2 / s^2 (CODATA 2010 value)
    const double electron_volt=1.602176565e-12;
    /// Electron mass in g (CODATA 2010 value)
    const double mass_electron=9.10938291e-28;
    /// Proton mass in g (CODATA 2010 value)
    const double mass_proton=1.672621777e-24;
    /// Neutron mass in g (CODATA 2010 value)
    const double mass_neutron=1.674927351e-24;
    /// Alpha particle mass in kg (CODATA 2010 value)
    const double mass_alpha=6.64465675e-24;
    /// Boltzmann constant in g cm^2 / K s^2 (CODATA 2010 value)
    const double boltzmann=1.3806488e-16;
    /// Atomic mass constant in g (CODATA 2010 value)
    const double unified_atomic_mass=1.660538921e-24;
    /// g
    const double solar_mass=1.9884e33;
    /// cm
    const double bohr_radius=5.291772083e-9;
    /// g / K^4 s^3
    const double stefan_boltzmann_constant=5.67039934436e-5;
    /// cm^2
    const double thomson_cross_section=6.65245853542e-25;

};



















