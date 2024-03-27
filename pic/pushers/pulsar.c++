#include <cmath> 

#include "pic/pushers/pulsar.h"
#include "tools/signum.h"
#include "external/iter/iter.h"
#include "tools/lerp.h"

#ifdef GPU
#include <nvtx3/nvToolsExt.h> 
#endif

using toolbox::sign;
using toolbox::lerp;


// cubic linear interpolation of the staggered emf to the 
// location of ind+dx+dy+dz
inline auto interpolate_fields(
    const toolbox::Mesh<float_m, 3>& exM,
    const toolbox::Mesh<float_m, 3>& eyM,
    const toolbox::Mesh<float_m, 3>& ezM,
    const toolbox::Mesh<float_m, 3>& bxM,
    const toolbox::Mesh<float_m, 3>& byM,
    const toolbox::Mesh<float_m, 3>& bzM,
    size_t ind, double dx, double dy, double dz, // interpolation target location
    size_t iy, size_t iz                         // mesh sizes in y and z dir
    ) -> std::tuple<double, double, double, double, double, double>
{

  double c000, c100, c010, c110, c001, c101, c011, c111;

  //ex with Yee-grid staggering of 1,0,0
  c000 = 0.5*(exM(ind         ) + exM(ind-1      ));
  c100 = 0.5*(exM(ind         ) + exM(ind+1      ));
  c010 = 0.5*(exM(ind   +iy   ) + exM(ind-1+iy   ));
  c110 = 0.5*(exM(ind   +iy   ) + exM(ind+1+iy   ));
  c001 = 0.5*(exM(ind      +iz) + exM(ind-1   +iz));
  c101 = 0.5*(exM(ind      +iz) + exM(ind+1   +iz));
  c011 = 0.5*(exM(ind   +iy+iz) + exM(ind-1+iy+iz));
  c111 = 0.5*(exM(ind   +iy+iz) + exM(ind+1+iy+iz));
  auto ex = lerp(c000, c100, c010, c110, c001, c101, c011, c111, dx, dy, dz);

  //ey 0,1,0
  c000 = 0.5*(eyM(ind        ) + eyM(ind  -iy   ));
  c100 = 0.5*(eyM(ind+1      ) + eyM(ind+1-iy   ));
  c010 = 0.5*(eyM(ind        ) + eyM(ind  +iy   ));
  c110 = 0.5*(eyM(ind+1      ) + eyM(ind+1+iy   ));
  c001 = 0.5*(eyM(ind     +iz) + eyM(ind  -iy+iz));
  c101 = 0.5*(eyM(ind+1   +iz) + eyM(ind+1-iy+iz));
  c011 = 0.5*(eyM(ind     +iz) + eyM(ind  +iy+iz));
  c111 = 0.5*(eyM(ind+1   +iz) + eyM(ind+1+iy+iz));
  auto ey = lerp(c000, c100, c010, c110, c001, c101, c011, c111, dx, dy, dz);

  //ez 0,0,1
  c000 = 0.5*(ezM(ind        ) + ezM(ind     -iz));
  c100 = 0.5*(ezM(ind+1      ) + ezM(ind+1   -iz));
  c010 = 0.5*(ezM(ind  +iy   ) + ezM(ind  +iy-iz));
  c110 = 0.5*(ezM(ind+1+iy   ) + ezM(ind+1+iy-iz));
  c001 = 0.5*(ezM(ind        ) + ezM(ind     +iz));
  c101 = 0.5*(ezM(ind+1      ) + ezM(ind+1   +iz));
  c011 = 0.5*(ezM(ind  +iy   ) + ezM(ind  +iy+iz));
  c111 = 0.5*(ezM(ind+1+iy   ) + ezM(ind+1+iy+iz));
  auto ez = lerp(c000, c100, c010, c110, c001, c101, c011, c111, dx, dy, dz);

  //-------------------------------------------------- 
  // bx 0,1,1
  c000 = 0.25*(bxM(ind  ) + bxM(ind-iy  ) + bxM(ind     -iz) +bxM(ind-iy-iz));
  c100 = 0.25*(bxM(ind+1) + bxM(ind+1-iy) + bxM(ind+1   -iz) +bxM(ind+1-iy-iz));
  c001 = 0.25*(bxM(ind  ) + bxM(ind+iz  ) + bxM(ind  -iy   ) +bxM(ind-iy+iz));
  c101 = 0.25*(bxM(ind+1) + bxM(ind+1+iz) + bxM(ind+1-iy   ) +bxM(ind+1-iy+iz));
  c010 = 0.25*(bxM(ind  ) + bxM(ind+iy  ) + bxM(ind     -iz) +bxM(ind+iy-iz));
  c110 = 0.25*(bxM(ind+1) + bxM(ind+1-iz) + bxM(ind+1+iy-iz) +bxM(ind+1+iy));
  c011 = 0.25*(bxM(ind  ) + bxM(ind+iy  ) + bxM(ind  +iy+iz) +bxM(ind+iz));
  c111 = 0.25*(bxM(ind+1) + bxM(ind+1+iy) + bxM(ind+1+iy+iz) +bxM(ind+1+iz));
  auto bx = lerp(c000, c100, c010, c110, c001, c101, c011, c111, dx, dy, dz);

  // by 1,0,1
  c000 = 0.25*(byM(ind-1-iz)+    byM(ind-1)+       byM(ind-iz)+      byM(ind));
  c100 = 0.25*(byM(ind-iz)+      byM(ind)+         byM(ind+1-iz)+    byM(ind+1));
  c001 = 0.25*(byM(ind-1)+       byM(ind-1+iz)+    byM(ind)+         byM(ind+iz));
  c101 = 0.25*(byM(ind)+         byM(ind+iz)+      byM(ind+1)+       byM(ind+1+iz));
  c010 = 0.25*(byM(ind-1+iy-iz)+ byM(ind-1+iy)+    byM(ind+iy-iz)+   byM(ind+iy));
  c110 = 0.25*(byM(ind+iy-iz)+   byM(ind+iy)+      byM(ind+1+iy-iz)+ byM(ind+1+iy));
  c011 = 0.25*(byM(ind-1+iy)+    byM(ind-1+iy+iz)+ byM(ind+iy)+      byM(ind+iy+iz));
  c111 = 0.25*(byM(ind+iy)+      byM(ind+iy+iz)+   byM(ind+1+iy)+    byM(ind+1+iy+iz));
  auto by = lerp(c000, c100, c010, c110, c001, c101, c011, c111, dx, dy, dz);

  // bz 1,1,0
  c000 = 0.25*(bzM(ind-1-iy)+    bzM(ind-1)+       bzM(ind-iy)+      bzM(ind));
  c100 = 0.25*(bzM(ind-iy)+      bzM(ind)+         bzM(ind+1-iy)+    bzM(ind+1));
  c001 = 0.25*(bzM(ind-1-iy+iz)+ bzM(ind-1+iz)+    bzM(ind-iy+iz)+   bzM(ind+iz));
  c101 = 0.25*(bzM(ind-iy+iz)+   bzM(ind+iz)+      bzM(ind+1-iy+iz)+ bzM(ind+1+iz));
  c010 = 0.25*(bzM(ind-1)+       bzM(ind-1+iy)+    bzM(ind)+         bzM(ind+iy));
  c110 = 0.25*(bzM(ind)+         bzM(ind+iy)+      bzM(ind+1)+       bzM(ind+1+iy));
  c011 = 0.25*(bzM(ind-1+iz)+    bzM(ind-1+iy+iz)+ bzM(ind+iz)+      bzM(ind+iy+iz));
  c111 = 0.25*(bzM(ind+iz)+      bzM(ind+iy+iz)+   bzM(ind+1+iz)+    bzM(ind+1+iy+iz));
  auto bz = lerp(c000, c100, c010, c110, c001, c101, c011, c111, dx, dy, dz);

  return {ex, ey, ez, bx, by, bz};
}



//ExB in units of c with E.B != 0 correction
//
// Based on Beklemieshev & Tessarotto 1999 and assumes B > E
inline auto ExB_drift_rel_approx( 
            double  ex,  double  ey,  double  ez,
            double  bx,  double  by,  double  bz
                     ) -> std::tuple<double, double, double, double, double>
{
  const double b2 = bx*bx + by*by + bz*bz;
  const double e2 = ex*ex + ey*ey + ez*ez;

  double vex = (ey*bz - ez*by)/(b2 + e2 + EPS);
  double vey = (ez*bx - ex*bz)/(b2 + e2 + EPS);
  double vez = (ex*by - ey*bx)/(b2 + e2 + EPS);
  double we2  = vex*vex + vey*vey + vez*vez; //|u|^2
  we2 = std::min(0.25, we2); // prevent NaN/overflow

  //// we -> ve
  double ginv = (1.0 - sqrt(1.0 - 4.0*we2 + EPS))/(2.0*we2 + EPS);
  vex *= ginv;
  vey *= ginv;
  vez *= ginv;

  double ve2 = vex*vex + vey*vey + vez*vez; //|v|^2
  double kappa = 1.0/(sqrt(1.0 - ve2) + EPS); // gamma factor

  return {vex, vey, vez, kappa, we2};
}


// b*: E.B corrected "relativistic" unit B field vector
inline auto mag_unit_vec_rel_approx( 
            double  ex,  double  ey,  double  ez,
            double  bx,  double  by,  double  bz,
            double we2
                     ) -> std::tuple< double, double, double>
{
  const double b2 = bx*bx + by*by + bz*bz;
  const double e2 = ex*ex + ey*ey + ez*ez;

  double edotb = ex*bx + ey*by + ez*bz;
  double eperpx = ex - edotb*bx/(b2 + EPS);
  double eperpy = ey - edotb*by/(b2 + EPS);
  double eperpz = ez - edotb*bz/(b2 + EPS);
  double eperp2 = eperpx*eperpx + eperpy*eperpy + eperpz*eperpz;


  double bp2 = 0.5*(b2 - e2 + (e2 + b2)*sqrt(1.0 - 4.0*we2)); // eq6
  double ep = edotb/(sqrt(bp2) + EPS); //eq 5

  double psi = ( ep/sqrt(bp2) )*(b2 - bp2)/(eperp2 + EPS);
  double eta = 1.0/(sqrt(b2)*sqrt( psi*psi*eperp2/(b2 + EPS) + 1.0) + EPS);
  double zeta = psi*eta;

  // rotation giving b* 
  double bnx = zeta*eperpx + eta*bx;
  double bny = zeta*eperpy + eta*by;
  double bnz = zeta*eperpz + eta*bz;

  // normalize to unit vector
  double bn = sqrt( bnx*bnx + bny*bny + bnz*bnz );
  bnx *= 1.0/(bn + EPS);
  bny *= 1.0/(bn + EPS);
  bnz *= 1.0/(bn + EPS);

  return {bnx, bny, bnz};
}


template<size_t D, size_t V>
void pic::PulsarPusher<D,V>::push_container(
    pic::ParticleContainer<D>& con, 
    pic::Tile<D>& tile
    )
{

#ifdef GPU
  nvtxRangePush(__PRETTY_FUNCTION__);
#endif

  const double c  = tile.cfl;
  const double qm = sign(con.q)/con.m; // q_s/m_s (sign only because emf are in units of q)
  const double m  = con.m; //mass

  // emf at the grid
  auto& gs = tile.get_grids(); 
  auto& exM = gs.ex;
  auto& eyM = gs.ey;
  auto& ezM = gs.ez;

  auto& bxM = gs.bx;
  auto& byM = gs.by;
  auto& bzM = gs.bz;

  //const int Nx = tile.mesh_lengths[0];
  //const int Ny = tile.mesh_lengths[1];
  //const int Nz = tile.mesh_lengths[2];

  // mesh sizes for 1D indexing
  const size_t iy = D >= 2 ? gs.ex.indx(0,1,0) - gs.ex.indx(0,0,0) : 0;
  const size_t iz = D >= 3 ? gs.ex.indx(0,0,1) - gs.ex.indx(0,0,0) : 0;

  auto mins = tile.mins;
  //auto maxs = tile.maxs;

  double ii=0.0, jj=0.0, kk=0.0;
  double dx=0.0, dy=0.0, dz=0.0;

  // loop over particles
  //UniIter::iterate([=] DEVCALLABLE (size_t n, pic::ParticleContainer<D>& con){
  for(size_t n=0; n<con.size(); n++){

    // particle location (assuming location coincides with gyro-center location)
    const auto R0x = con.loc(0,n);
    const auto R0y = con.loc(1,n);
    const auto R0z = con.loc(2,n);

    // particle velocity
    const auto vel0n = con.vel(0,n);
    const auto vel1n = con.vel(1,n);
    const auto vel2n = con.vel(2,n);

    //--------------------------------------------------
    // cell integer locations
    if(D > 0) ii = floor(R0x);
    if(D > 1) jj = floor(R0y);
    if(D > 2) kk = floor(R0z);

    // cell shifts
    if(D > 0) dx = R0x - ii;
    if(D > 1) dy = R0y - jj;
    if(D > 2) dz = R0z - kk;

    // normalize to tile units
    if(D > 0) ii -= mins[0];
    if(D > 1) jj -= mins[1];
    if(D > 2) kk -= mins[2];
      
    //--------------------------------------------------
    //internal field interpolation 

    // 1D reference index
    const size_t ind0 = gs.ex.indx(ii,jj,kk);

    auto [ex0, ey0, ez0, bx0, by0, bz0] = interpolate_fields(
                                      exM, eyM, ezM, bxM, byM, bzM, 
                                      ind0, dx, dy, dz, 
                                      iy, iz);

    // add external field component from pusher
    //TODO: why cinv here in E?
    ex0 += this->get_ex_ext(0,0,0); 
    ey0 += this->get_ey_ext(0,0,0);
    ez0 += this->get_ez_ext(0,0,0);
    bx0 += this->get_bx_ext(0,0,0);
    by0 += this->get_by_ext(0,0,0);
    bz0 += this->get_bz_ext(0,0,0);

    // normalize units to be ready for the internal scheme
    ex0 *= 1.0/c;
    ey0 *= 1.0/c;
    ez0 *= 1.0/c;    
    bx0 *= 1.0/c;    
    by0 *= 1.0/c;
    bz0 *= 1.0/c;

    // field strenghts
    const double b0 = sqrt( bx0*bx0 + by0*by0 + bz0*bz0 );
    const double e0 = sqrt( ex0*ex0 + ey0*ey0 + ez0*ez0 );

    //-------------------------------------------------- 
    // relativistic ExB drift velocity
    auto [vex0, vey0, vez0, kappa0, we2] = ExB_drift_rel_approx( ex0, ey0, ez0, bx0, by0, bz0 );

    // magnetic field unit vector b
    auto [bnx0, bny0, bnz0] = mag_unit_vec_rel_approx( ex0, ey0, ez0, bx0, by0, bz0, we2);

    // epar = e.b
    double epar = ex0*bnx0 + ey0*bny0 + ez0*bnz0;

    // project velocity to the b field; upar at t_n-1/2 via u_par = (u . b)b
    double upar01  = vel0n*bnx0 + vel1n*bny0 + vel2n*bnz0;

    //--------------------------------------------------
    // gyro motion four-velocity

    // NOTE: assuming full gamma factor here from previous step
    double G0 = sqrt(1.0 + vel0n*vel0n + vel1n*vel1n + vel2n*vel2n );

    double ugx = vel0n - upar01*bnx0 - vex0*G0;
    double ugy = vel1n - upar01*bny0 - vey0*G0;
    double ugz = vel2n - upar01*bnz0 - vez0*G0;
    double ug2 = ugx*ugx + ugy*ugy + ugz*ugz;

    // magnetic moment = m u_g^2/2 B_0 \gamma
    double mu = m*ug2/(2.0*b0*kappa0);

    //--------------------------------------------------
    // upar at t_n+1/2
    // increase velocity with parallel E field: u += q/m dt Epar
    // NOTE: standard c * vel + qm*epar changed to vel + qm*epar/c
    // NOTE: cinv is multiplied to b0 in the beginning
    upar01 += qm*epar;
    const double k0 = sqrt(1.0 + upar01*upar01 + ug2 ); // gamma


    //--------------------------------------------------
    // preparing for step1; assume that solution at t=1 is at the same cell 
    double R1x = R0x;
    double R1y = R0y;
    double R1z = R0z;

    //-------------------------------------------------- 
    // implicit iteration (done in N steps

    double vn1x, vn1y, vn1z, un1x, un1y, un1z, G1;

    for(size_t iter=0; iter<3; iter++){

      if(D > 0) ii = floor(R1x);
      if(D > 1) jj = floor(R1y);
      if(D > 2) kk = floor(R1z);

      // cell shifts
      if(D > 0) dx = R1x - ii;
      if(D > 1) dy = R1y - jj;
      if(D > 2) dz = R1z - kk;

      // normalize to tile units
      if(D > 0) ii -= mins[0];
      if(D > 1) jj -= mins[1];
      if(D > 2) kk -= mins[2];
      
      const size_t ind1 = gs.ex.indx(ii,jj,kk); // 1D reference index

      auto [ex1, ey1, ez1, bx1, by1, bz1] = interpolate_fields(
                                      exM, eyM, ezM, bxM, byM, bzM, 
                                      ind1, dx, dy, dz, 
                                      iy, iz);

      // add external field component from pusher
      //TODO: why cinv here in E?
      ex1 += this->get_ex_ext(0,0,0); 
      ey1 += this->get_ey_ext(0,0,0);
      ez1 += this->get_ez_ext(0,0,0);
      bx1 += this->get_bx_ext(0,0,0);
      by1 += this->get_by_ext(0,0,0);
      bz1 += this->get_bz_ext(0,0,0);

      // normalize units to be ready for the internal scheme
      ex1 *= 1.0/c;
      ey1 *= 1.0/c;
      ez1 *= 1.0/c;    
      bx1 *= 1.0/c;    
      by1 *= 1.0/c;
      bz1 *= 1.0/c;

      // field strenghts
      const double b1 = sqrt( bx1*bx1 + by1*by1 + bz1*bz1 );
      //const double e1 = sqrt( ex1*ex1 + ey1*ey1 + ez1*ez1 );
    
      //-------------------------------------------------- 
      //ExB in units of c at new location
      auto [vex1, vey1, vez1, kappa1, we2] = ExB_drift_rel_approx( ex1, ey1, ez1, bx1, by1, bz1 );

      //-------------------------------------------------- 
      // magnetic field unit vector b at new location
      auto [bnx1, bny1, bnz1] = mag_unit_vec_rel_approx( ex1, ey1, ez1, bx1, by1, bz1, we2);

      //-------------------------------------------------- 
      // location update
      // NOTE: assuming mu -> 0 zero during the time step

      // Construct new gyrovelocity from conservation of magnetic moment
      // magnetic moment = m u_g^2/2 B_0 \gamma
        
      mu = 0.0; // NOTE: synchrotron losses are assumed to bring mag. mom. to zero
      double ug2n = mu*2.0*b1*kappa1/m;

      const double k1 = sqrt(1.0 + upar01*upar01 + ug2n);     // gamma

      // update Gamma estimates
      G0 = k0*kappa0; // inv Gamma at t = n FIXME: or previous G0?
      G1 = k1*kappa1; // inv Gamma at t = n+1

      // GCA coordinate-velocity at v_n+1/2 
      vn1x = 0.5*upar01*( bnx0/G0 + bnx1/G1 ) + 0.5*(vex0 + vex1);
      vn1y = 0.5*upar01*( bny0/G0 + bny1/G1 ) + 0.5*(vey0 + vey1);
      vn1z = 0.5*upar01*( bnz0/G0 + bnz1/G1 ) + 0.5*(vez0 + vez1);

      //-------------------------------------------------- 
      // newest GCA four velocity at u_n+1 (for velocity update)
      // gyromotion is taken to be same direction as in the previous step (hence ug_i/sqrt(ug2))
      un1x = upar01*bnx1 + vex1*G1 + ug2n*ugx/sqrt(ug2);
      un1y = upar01*bny1 + vey1*G1 + ug2n*ugy/sqrt(ug2);
      un1z = upar01*bnz1 + vez1*G1 + ug2n*ugz/sqrt(ug2);

      //-------------------------------------------------- 
      // location error
      double Hx = R1x - (R0x + c*vn1x);
      double Hy = R1y - (R0y + c*vn1y);
      double Hz = R1z - (R0z + c*vn1z);
      double H = sqrt(Hx*Hx + Hy*Hy + Hz*Hz);

      //-------------------------------------------------- 
      // guiding center location update
      if(D>=1) R1x = R0x + vn1x*c;
      if(D>=2) R1y = R0y + vn1y*c;
      if(D>=3) R1z = R0z + vn1z*c;


#ifdef DEBUG
      //--------------------------------------------------
      // debug
      bool debug_flag2 = 
        std::isnan(un1x) ||
        std::isnan(un1y) ||
        std::isnan(un1z) ||
        std::isnan(R1x) ||
        std::isnan(R1y) ||
        std::isnan(R1z);   

      if(debug_flag2) {
          
        double b1 = sqrt( bx1*bx1 + by1*by1 + bz1*bz1 );
        double e1 = sqrt( ex1*ex1 + ey1*ey1 + ez1*ez1 );

        // E.B violation; i.e .how bad is the pure ExB drift assumption
        double EB0 = (ex0*bx0 + ey0*by0 + ez0*bz0)/b0/b0;
        double EB1 = (ex1*bx1 + ey1*by1 + ez1*bz1)/b1/b1;
        double bn0 = bnx0*bnx0 + bny0*bny0 + bnz0*bnz0;
        double bn1 = bnx1*bnx1 + bny1*bny1 + bnz1*bnz1;

        double egtb0 = e0/b0;
        double egtb1 = e1/b1;

        std::cout 
        << "iter:" << iter << " H:" << H << " E.B:" << EB0 << " " << EB1
        << " e/b0:" << egtb0 << " e/b1:" << egtb1
        << "\n"
        << " R0x:" << R0x << " R0y:" << R0y << " R0z:" << R0z 
        << "\n"
        << " R1x:" << R1x << " R1y:" << R1y << " R1z:" << R1z  
        << "\n"
        << " Hx:"  << Hx  << " Hy:"  << Hy  << " Hz:"  << Hz
        << "\n"
        << " ex0:" << ex0 << " ey0:" << ey0 << " ez0:" << ez0
        << " bx0:" << bx0 << " by0:" << by0 << " bz0:" << bz0
        << " bnx0:" << bnx0 << " bny:" << bny0 << " bnz:" << bnz0 << " b: " << bn0
        //<< "\n"
        //<< " *bnx0:" << bnx0_2 << " bny:" << bny0_2 << " bnz:" << bnz0_2
        << "\n"
        << " bnx1:" << bnx1 << " bny:" << bny1 << " bnz:" << bnz1 << " b: " << bn1
        << "\n"
        //<< " eparx:" << eparx << " epary:" << epary << " eparz:" << eparz
        //<< " uparx:" << uparx01 << " upary:" << upary01 << " uparz:" << uparz01 << 
        " upar:" << upar01
        << " k0:" << k0 
        << " k1:" << k1 
        << " G0:" << G0 
        << " G1:" << G1 
        << "\n"
        << " mu:" << mu 
        << " ug:" << sqrt(ug2) << " ugn:" << sqrt(ug2n)
        << "\n"
        << " vex0:" << vex0 << " vey:" << vey0 << " vez:" << vez0 << " kappa0:" << kappa0
        << "\n"
        << " vex1:" << vex1 << " vey:" << vey1 << " vez:" << vez1 << " kappa1:" << kappa1
        << "\n"
        << " i:" << ii << " j:" << jj << " k:" << kk
        << " dx:" << dx << " dy:" << dy << " dz:" << dz
        << "\n\n";
        std::cout << std::flush;

        assert(false);
      }
#endif

    } // end of loop over particles
    //}, con.size(), con);

    // converged

    // update velocity
    con.vel(0,n) = un1x;
    con.vel(1,n) = un1y;
    con.vel(2,n) = un1z;

    // additional gravity

    // location at the half time step
    double xx = 0.5*(R0x + R1x); 
    double yy = 0.5*(R0y + R1y); 
    double zz = 0.5*(R0z + R1z); 

    // velocity at the half time step 
    //double uxt = 0.5*(un1x + vel0n);
    //double uyt = 0.5*(un1y + vel1n);
    //double uzt = 0.5*(un1z + vel2n);
    //double gamt = sqrt(1.0 + uxt*uxt + uyt*uyt + uzt*uzt); // gamma at half time step

    // distance from the star 
    double rad = sqrt( pow(xx-cenx, 2) + pow(yy-ceny, 2) + pow(zz-cenz, 2) );

    // gravitational force is F_g = G M m/r^2
    // surface gravity is given as g = G M /r^2 
    // hence gravitational acceleration is a_grav = F_grav/m = G M/r^2 = g_0 (R/r)^2
    double a_grav = gravity_const*pow(rad_star/rad, 2);

    // gravitational acceleration to the radial direction (towards star)  
    double gravx = a_grav*( cenx - xx )/rad;
    double gravy = a_grav*( ceny - yy )/rad;
    double gravz = a_grav*( cenz - zz )/rad;

    // TODO c = dt here or not?
    con.vel(0,n) += c*gravx;
    con.vel(1,n) += c*gravy;
    con.vel(2,n) += c*gravz;

    // position update from iteration, new location is the gyro center position
    if(D>=1) con.loc(0,n) = R1x;
    if(D>=2) con.loc(1,n) = R1y;
    if(D>=3) con.loc(2,n) = R1z;  


#ifdef DEBUG
    //std::cout << "PulsarPusher: g0" << gravity_const << " a_g: " << a_grav << " gravy " << gravy << " \n";
      
    //-------------------------------------------------- 
    // debug

    bool debug_flag = 
    std::isnan(con.vel(0,n)) ||
    std::isnan(con.vel(1,n)) ||
    std::isnan(con.vel(2,n)) ||
    std::isnan(con.loc(0,n)) ||
    std::isnan(con.loc(1,n)) ||
    std::isnan(con.loc(2,n));   

    // test that we do not move too far
    double dxp = std::abs( con.loc(0,n) - R1x );
    double dyp = std::abs( con.loc(1,n) - R1y );
    double dzp = std::abs( con.loc(2,n) - R1z );
    if( (dxp > c) || (dyp > c) || (dzp > c) ) debug_flag = true;

    //if(1./kinv01 > 30.0) debug_flag = true;
    //if(true){
    if(debug_flag){
      std::cout 
        << " n:" << n
        << " loc0n:" << con.loc(0,n) << " loc1n:" << con.loc(1,n) << " loc2n:" << con.loc(2,n)
        << " dx:" << dxp << " " << dyp << " " << dzp 
        << " ex0:" << ex0 << " ey0:" << ey0 << " ez0:" << ez0
        << " bx0:" << bx0 << " by0:" << by0 << " bz0:" << bz0
        << " bnx0:" << bnx0 << " bny:" << bny0 << " bnz:" << bnz0
        //<< " bnx1:" << bnx1 << " bny:" << bny1 << " bnz:" << bnz1
        //<< " eparx:" << eparx << " epary:" << epary << " eparz:" << eparz
        //<< " uparx:" << uparx01 << " upary:" << upary01 << " uparz:" << uparz01
        << " upar:" << upar01
        << " k0:" << k0 
        //<< " k1:" << k1 
        << " vex0:" << vex0 << " vey:" << vey0 << " vez:" << vez0 << " kappa0:" << kappa0
        //<< " vex1:" << vex1 << " vey:" << vey1 << " vez:" << vez1 << " kappa1:" << kappa1
        << "\n";
      std::cout << std::flush;

      assert(false);
    }
#endif
  //--------------------------------------------------


  } // end of loop over containers


  UniIter::sync();

#ifdef GPU
  nvtxRangePop();
#endif
}



//--------------------------------------------------
// explicit template instantiation

template class pic::PulsarPusher<1,3>; // 1D3V
template class pic::PulsarPusher<2,3>; // 2D3V
template class pic::PulsarPusher<3,3>; // 3D3V

