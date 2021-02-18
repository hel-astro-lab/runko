#include "rgca.h"

#include <cmath> 
#include "../../tools/signum.h"

using toolbox::sign;



inline real_long _lerp(
      real_long c000,
      real_long c100,
      real_long c010,
      real_long c110,
      real_long c001,
      real_long c101,
      real_long c011,
      real_long c111,
      real_long dx, real_long dy, real_long dz
      ) 
{
      real_long c00 = c000 * (1.0-dx) + c100 * dx;
      real_long c10 = c010 * (1.0-dx) + c110 * dx;
      real_long c0  = c00  * (1.0-dy) + c10  * dy;
      real_long c01 = c001 * (1.0-dx) + c101 * dx;
      real_long c11 = c011 * (1.0-dx) + c111 * dx;
      real_long c1  = c01  * (1.0-dy) + c11  * dy;
      real_long c   = c0   * (1.0-dz) + c1   * dz;
      return c;
}


//-------------------------------------------------- 
inline auto ExB_drift( 
            real_long  ex,  real_long  ey,  real_long  ez,
            real_long  bx,  real_long  by,  real_long  bz
                     ) -> std::tuple<real_long, real_long, real_long, real_long, real_long>
{
    const real_long b = sqrt( bx*bx + by*by + bz*bz );

    real_long vex = (ey*bz - ez*by)/(b*b + EPS);
    real_long vey = (ez*bx - ex*bz)/(b*b + EPS);
    real_long vez = (ex*by - ey*bx)/(b*b + EPS);

    real_long ve2 = vex*vex + vey*vey + vez*vez; //|u|^2
    real_long kappa = 1.0/(sqrt(1. - ve2) + EPS); // gamma factor

    return {vex, vey, vez, kappa, 0.};
}
    

//ExB in units of c with E.B !=  correction
inline auto ExB_drift_rel( 
            real_long  ex,  real_long  ey,  real_long  ez,
            real_long  bx,  real_long  by,  real_long  bz
                     ) -> std::tuple<real_long, real_long, real_long, real_long, real_long>
{
    const real_long b = sqrt( bx*bx + by*by + bz*bz );
    const real_long e = sqrt( ex*ex + ey*ey + ez*ez );

    real_long vex = (ey*bz - ez*by)/(b*b + e*e + EPS);
    real_long vey = (ez*bx - ex*bz)/(b*b + e*e + EPS);
    real_long vez = (ex*by - ey*bx)/(b*b + e*e + EPS);
    real_long we2  = vex*vex + vey*vey + vez*vez; //|u|^2

    //// we -> ve
    real_long ginv = (1. - sqrt(1. - 4.*we2))/(2.*we2 + EPS);
    vex *= ginv;
    vey *= ginv;
    vez *= ginv;

    real_long ve2 = vex*vex + vey*vey + vez*vez; //|u|^2
    real_long kappa = 1.0/(sqrt(1. - ve2) + EPS); // gamma factor

    return {vex, vey, vez, kappa, we2};
}


// b: normal unit vector 
inline auto mag_unit_vec( 
            real_long  bx,  real_long  by,  real_long  bz
                     ) -> std::tuple< real_long, real_long, real_long>
{
    const real_long b = sqrt( bx*bx + by*by + bz*bz );

    real_long bnx = bx/(b + EPS); 
    real_long bny = by/(b + EPS); 
    real_long bnz = bz/(b + EPS); 

    return {bnx, bny, bnz};
}


// b*: exact, E.B corrected "relativistic" unit B field vector
inline auto mag_unit_vec_rel( 
            real_long  ex,  real_long  ey,  real_long  ez,
            real_long  bx,  real_long  by,  real_long  bz,
            real_long we2
                     ) -> std::tuple< real_long, real_long, real_long>
{
    const real_long b = sqrt( bx*bx + by*by + bz*bz );
    const real_long e = sqrt( ex*ex + ey*ey + ez*ez );

    real_long edotb = ex*bx + ey*by + ez*bz;
    real_long eperpx = ex - edotb*bx/(b*b + EPS);
    real_long eperpy = ey - edotb*by/(b*b + EPS);
    real_long eperpz = ez - edotb*bz/(b*b + EPS);
    real_long eperp = eperpx*eperpx + eperpy*eperpy + eperpz*eperpz;

    real_long bp = sqrt(0.5*(b*b - e*e + (e*e + b*b)*sqrt(1.0 - 4.*we2))); // eq6
    real_long ep = edotb/(bp + EPS); //eq 5

    real_long psi = ep*(b*b - bp*bp)/(bp*eperp*eperp + EPS);
    real_long eta = 1.0/(b*sqrt( psi*psi*eperp*eperp/(b*b + EPS) + 1.) + EPS);
    real_long zeta = psi*eta;

    // rotation giving b* 
    real_long bnx = zeta*eperpx + eta*bx;
    real_long bny = zeta*eperpy + eta*by;
    real_long bnz = zeta*eperpz + eta*bz;

    // normalize to unit vector
    real_long bn = sqrt( bnx*bnx + bny*bny + bnz*bnz );
    bnx *= 1.0/(bn + EPS);
    bny *= 1.0/(bn + EPS);
    bnz *= 1.0/(bn + EPS);

    return {bnx, bny, bnz};
}



template<size_t D, size_t V>
void pic::rGCAPusher<D,V>::push_container(
    pic::ParticleContainer<D>& container, 
    pic::Tile<D>& tile
    )
{
  int nparts = container.size();

  // initialize pointers to particle arrays
  real_prtcl* loc[3];
  for( int i=0; i<3; i++) loc[i] = &( container.loc(i,0) );

  real_prtcl* vel[3];
  for( int i=0; i<3; i++) vel[i] = &( container.vel(i,0) );


  real_long ex0 = 0.0, ey0 = 0.0, ez0 = 0.0;
  real_long bx0 = 0.0, by0 = 0.0, bz0 = 0.0;

  // make sure E and B tmp arrays are of correct size
  if(container.Epart.size() != (size_t)3*nparts)
    container.Epart.resize(3*nparts);
  if(container.Bpart.size() != (size_t)3*nparts)
    container.Bpart.resize(3*nparts);

  // fields at prtcl loc
  real_prtcl *exP, *eyP, *ezP, *bxP, *byP, *bzP;
  exP = &( container.Epart[0*nparts] );
  eyP = &( container.Epart[1*nparts] );
  ezP = &( container.Epart[2*nparts] );

  bxP = &( container.Bpart[0*nparts] );
  byP = &( container.Bpart[1*nparts] );
  bzP = &( container.Bpart[2*nparts] );

  const int Nx = tile.mesh_lengths[0];
  const int Ny = tile.mesh_lengths[1];
  const int Nz = tile.mesh_lengths[2];

  // fields at grid
  auto& yee = tile.get_yee();

  auto& exM = yee.ex;
  auto& eyM = yee.ey;
  auto& ezM = yee.ez;

  auto& bxM = yee.bx;
  auto& byM = yee.by;
  auto& bzM = yee.bz;

  // loop over particles
  int n1 = 0;
  int n2 = nparts;

  real_long c = tile.cfl;
  real_long cinv = 1.0/c; 

  // half charge-to-mass ratio (sign only because fields are in units of q)
  real_long qm = sign(container.q)/container.m;
  real_long me = container.m;

  real_long loc0n, loc1n, loc2n;
  real_long vel0n, vel1n, vel2n;

  //real_long edotb, eperpx, eperpy, eperpz, eperp, bp, ep, psi, eta, zeta;
  real_long ugx, ugy, ugz, ug2, ug2n;
  real_long G0, G1;
  real_long mu;

  // work variables
  real_long bx1, by1, bz1;
  real_long ex1, ey1, ez1;

  real_long vn1x, vn1y, vn1z;
  real_long un1x, un1y, un1z;

  real_long c000, c100, c010, c110, c001, c101, c011, c111;

  // mesh sizes for 1D indexing
  const size_t iy = D >= 2 ? yee.ex.indx(0,1,0) - yee.ex.indx(0,0,0) : 0;
  const size_t iz = D >= 3 ? yee.ex.indx(0,0,1) - yee.ex.indx(0,0,0) : 0;
  auto mins = tile.mins;

  real_long dx=0.0, dy=0.0, dz=0.0;
  int i=0, j=0, k=0;

  // loop over prtcls
  for(int n=n1; n<n2; n++) {
      bool crash_flag = false;

    loc0n = static_cast<real_long>( loc[0][n] );
    loc1n = static_cast<real_long>( loc[1][n] );
    loc2n = static_cast<real_long>( loc[2][n] );

    vel0n = static_cast<real_long>( vel[0][n] );
    vel1n = static_cast<real_long>( vel[1][n] );
    vel2n = static_cast<real_long>( vel[2][n] );

    // read particle-specific fields
    ex0 = static_cast<real_long>( (exP[n] + this->get_ex_ext(0,0,0))*cinv );
    ey0 = static_cast<real_long>( (eyP[n] + this->get_ey_ext(0,0,0))*cinv );
    ez0 = static_cast<real_long>( (ezP[n] + this->get_ez_ext(0,0,0))*cinv );
    bx0 = static_cast<real_long>( (bxP[n] + this->get_bx_ext(0,0,0))*cinv );
    by0 = static_cast<real_long>( (byP[n] + this->get_by_ext(0,0,0))*cinv );
    bz0 = static_cast<real_long>( (bzP[n] + this->get_bz_ext(0,0,0))*cinv );

    //-------------------------------------------------- 
    // iterate: step0
      
    const real_long R0x = loc0n;
    const real_long R0y = loc1n;
    const real_long R0z = loc2n;

    const real_long b0 = sqrt( bx0*bx0 + by0*by0 + bz0*bz0 );
    const real_long e0 = sqrt( ex0*ex0 + ey0*ey0 + ez0*ez0 );

    //-------------------------------------------------- 
    //ExB in units of c

    // non-rel / rel ExB drift velocity
    //auto [vex0, vey0, vez0, kappa0, we2] = ExB_drift( ex0, ey0, ez0, bx0, by0, bz0 );
    auto [vex0, vey0, vez0, kappa0, we2] = ExB_drift_rel( ex0, ey0, ez0, bx0, by0, bz0 );

    //-------------------------------------------------- 
    // magnetic field unit vector b
      
    //auto [bnx0, bny0, bnz0] = mag_unit_vec(bx0, by0, bz0);
    auto [bnx0, bny0, bnz0] = mag_unit_vec_rel( ex0, ey0, ez0, bx0, by0, bz0, we2);

    //--------------------------------------------------
    // epar = e.b
    real_long epar = ex0*bnx0 + ey0*bny0 + ez0*bnz0;

    // project velocity to the b field; upar at t_n-1/2 via u_par = (u . b)b
    real_long upar01  = vel0n*bnx0 + vel1n*bny0 + vel2n*bnz0;


    //--------------------------------------------------
    // gyro motion four-velocity
      
    // NOTE: assuming full gamma factor here from previous step
    G0 = sqrt(1.0 + vel0n*vel0n + vel1n*vel1n + vel2n*vel2n );

    ugx = vel0n - upar01*bnx0 - vex0*G0;
    ugy = vel1n - upar01*bny0 - vey0*G0;
    ugz = vel2n - upar01*bnz0 - vez0*G0;
    ug2 = ugx*ugx + ugy*ugy + ugz*ugz;

    // magnetic moment = m u_g^2/2 B_0 \gamma
    // FIXME: this is taken as kappa0 in paper eqs; correct?
    //mu = me*ug2/(2.0*b0*G);
    mu = me*ug2/(2.0*b0*kappa0);

    //--------------------------------------------------
    // upar at t_n+1/2
    // increase velocity with parallel E field: u += q/m dt Epar
    // NOTE: standard c * vel + qm*epar changed to vel + qm*epar/c
    // NOTE: cinv is multiplied to b0 in the beginning
    // FIXME: or multiply cinv here?
    upar01 += qm*epar;
    const real_long k0 = sqrt(1.0 + upar01*upar01 + ug2 );     // gamma

    //--------------------------------------------------
    // step1
    real_long R1x = R0x;
    real_long R1y = R0y;
    real_long R1z = R0z;


    for(size_t iter=0; iter<5; iter++){

      // re-use previous interpolation for first step
      if(iter == 0) {
        ex1 = ex0; ey1 = ey0; ez1 = ez0;
        bx1 = bx0; by1 = by0; bz1 = bz0;
      } else {

        if(D >= 1) i  = static_cast<int>(floor(R1x));
        if(D >= 2) j  = static_cast<int>(floor(R1y));
        if(D >= 3) k  = static_cast<int>(floor(R1z));

        if(D >= 1) dx = R1x - i;
        if(D >= 2) dy = R1y - j;
        if(D >= 3) dz = R1z - k;

        // normalize to tile units
        if(D >= 1) i -= mins[0];
        if(D >= 2) j -= mins[1];
        if(D >= 3) k -= mins[2];

        //extrapolate if on tile boundary
        // min left side of tile
        if(D >= 1) { if(i <= -2 ) { dx -= i+2; i = -2; } }
        if(D >= 2) { if(j <= -2 ) { dy -= j+2; j = -2; } }
        if(D >= 3) { if(k <= -2 ) { dz -= k+2; k = -2; } }
        
        // max right side of tile
        if(D >= 1) { if(i >= Nx+1 ) { dx += i-Nx-1; i = Nx+1; } }
        if(D >= 2) { if(j >= Ny+1 ) { dy += j-Ny-1; j = Ny+1; } }
        if(D >= 3) { if(k >= Nz+1 ) { dz += k-Nz-1; k = Nz+1; } }


        const size_t ind = yee.ex.indx(i,j,k);

        //ex
        c000 = 0.5*(exM(ind       ) +exM(ind-1      ));
        c100 = 0.5*(exM(ind       ) +exM(ind+1      ));
        c010 = 0.5*(exM(ind+iy    ) +exM(ind-1+iy   ));
        c110 = 0.5*(exM(ind+iy    ) +exM(ind+1+iy   ));
        c001 = 0.5*(exM(ind+iz    ) +exM(ind-1+iz   ));
        c101 = 0.5*(exM(ind+iz    ) +exM(ind+1+iz   ));
        c011 = 0.5*(exM(ind+iy+iz ) +exM(ind-1+iy+iz));
        c111 = 0.5*(exM(ind+iy+iz ) +exM(ind+1+iy+iz));
        ex1 = _lerp(c000, c100, c010, c110, c001, c101, c011, c111, dx, dy, dz);
        ex1 += this->get_ex_ext(0,0,0);

        //ey
        c000 = 0.5*(eyM(ind      ) +eyM(ind-iy     ));
        c100 = 0.5*(eyM(ind+1    ) +eyM(ind+1-iy   ));
        c010 = 0.5*(eyM(ind      ) +eyM(ind+iy     ));
        c110 = 0.5*(eyM(ind+1    ) +eyM(ind+1+iy   ));
        c001 = 0.5*(eyM(ind+iz   ) +eyM(ind-iy+iz  ));
        c101 = 0.5*(eyM(ind+1+iz ) +eyM(ind+1-iy+iz));
        c011 = 0.5*(eyM(ind+iz   ) +eyM(ind+iy+iz  ));
        c111 = 0.5*(eyM(ind+1+iz ) +eyM(ind+1+iy+iz));
        ey1 = _lerp(c000, c100, c010, c110, c001, c101, c011, c111, dx, dy, dz);
        ey1 += this->get_ey_ext(0,0,0);

        //ez
        c000 = 0.5*(ezM(ind      ) + ezM(ind-iz     ));
        c100 = 0.5*(ezM(ind+1    ) + ezM(ind+1-iz   ));
        c010 = 0.5*(ezM(ind+iy   ) + ezM(ind+iy-iz  ));
        c110 = 0.5*(ezM(ind+1+iy ) + ezM(ind+1+iy-iz));
        c001 = 0.5*(ezM(ind      ) + ezM(ind+iz     ));
        c101 = 0.5*(ezM(ind+1    ) + ezM(ind+1+iz   ));
        c011 = 0.5*(ezM(ind+iy   ) + ezM(ind+iy+iz  ));
        c111 = 0.5*(ezM(ind+1+iy ) + ezM(ind+1+iy+iz));
        ez1 = _lerp(c000, c100, c010, c110, c001, c101, c011, c111, dx, dy, dz);
        ez1 += this->get_ez_ext(0,0,0);

        //-------------------------------------------------- 
        // bx
        c000 = 0.25*( bxM(ind)+   bxM(ind-iy)+   bxM(ind-iz)+      bxM(ind-iy-iz));
        c100 = 0.25*( bxM(ind+1)+ bxM(ind+1-iy)+ bxM(ind+1-iz)+    bxM(ind+1-iy-iz));
        c001 = 0.25*( bxM(ind)+   bxM(ind+iz)+   bxM(ind-iy)+      bxM(ind-iy+iz));
        c101 = 0.25*( bxM(ind+1)+ bxM(ind+1+iz)+ bxM(ind+1-iy)+    bxM(ind+1-iy+iz));
        c010 = 0.25*( bxM(ind)+   bxM(ind+iy)+   bxM(ind-iz)+      bxM(ind+iy-iz));
        c110 = 0.25*( bxM(ind+1)+ bxM(ind+1-iz)+ bxM(ind+1+iy-iz)+ bxM(ind+1+iy));
        c011 = 0.25*( bxM(ind)+   bxM(ind+iy)+   bxM(ind+iy+iz)+   bxM(ind+iz));
        c111 = 0.25*( bxM(ind+1)+ bxM(ind+1+iy)+ bxM(ind+1+iy+iz)+ bxM(ind+1+iz));
        bx1 = _lerp(c000, c100, c010, c110, c001, c101, c011, c111, dx, dy, dz);
        bx1 += this->get_bx_ext(0,0,0);

        // by
        c000 = 0.25*( byM(ind-1-iz)+    byM(ind-1)+       byM(ind-iz)+      byM(ind));
        c100 = 0.25*( byM(ind-iz)+      byM(ind)+         byM(ind+1-iz)+    byM(ind+1));
        c001 = 0.25*( byM(ind-1)+       byM(ind-1+iz)+    byM(ind)+         byM(ind+iz));
        c101 = 0.25*( byM(ind)+         byM(ind+iz)+      byM(ind+1)+       byM(ind+1+iz));
        c010 = 0.25*( byM(ind-1+iy-iz)+ byM(ind-1+iy)+    byM(ind+iy-iz)+   byM(ind+iy));
        c110 = 0.25*( byM(ind+iy-iz)+   byM(ind+iy)+      byM(ind+1+iy-iz)+ byM(ind+1+iy));
        c011 = 0.25*( byM(ind-1+iy)+    byM(ind-1+iy+iz)+ byM(ind+iy)+      byM(ind+iy+iz));
        c111 = 0.25*( byM(ind+iy)+      byM(ind+iy+iz)+   byM(ind+1+iy)+    byM(ind+1+iy+iz));
        by1 = _lerp(c000, c100, c010, c110, c001, c101, c011, c111, dx, dy, dz);
        by1 += this->get_by_ext(0,0,0);

        // bz
        c000 = 0.25*( bzM(ind-1-iy)+    bzM(ind-1)+       bzM(ind-iy)+      bzM(ind));
        c100 = 0.25*( bzM(ind-iy)+      bzM(ind)+         bzM(ind+1-iy)+    bzM(ind+1));
        c001 = 0.25*( bzM(ind-1-iy+iz)+ bzM(ind-1+iz)+    bzM(ind-iy+iz)+   bzM(ind+iz));
        c101 = 0.25*( bzM(ind-iy+iz)+   bzM(ind+iz)+      bzM(ind+1-iy+iz)+ bzM(ind+1+iz));
        c010 = 0.25*( bzM(ind-1)+       bzM(ind-1+iy)+    bzM(ind)+         bzM(ind+iy));
        c110 = 0.25*( bzM(ind)+         bzM(ind+iy)+      bzM(ind+1)+       bzM(ind+1+iy));
        c011 = 0.25*( bzM(ind-1+iz)+    bzM(ind-1+iy+iz)+ bzM(ind+iz)+      bzM(ind+iy+iz));
        c111 = 0.25*( bzM(ind+iz)+      bzM(ind+iy+iz)+   bzM(ind+1+iz)+    bzM(ind+1+iy+iz));
        bz1 = _lerp(c000, c100, c010, c110, c001, c101, c011, c111, dx, dy, dz);
        bz1 += this->get_bz_ext(0,0,0);

        ex1 *= cinv;
        ey1 *= cinv;
        ez1 *= cinv;
        bx1 *= cinv;
        by1 *= cinv;
        bz1 *= cinv;
      }

      //-------------------------------------------------- 

      //-------------------------------------------------- 
      //ExB in units of c at new location

      // non-rel / rel ExB drift velocity at the new location
      //auto [vex1, vey1, vez1, kappa1, we2] = ExB_drift(     ex1, ey1, ez1, bx1, by1, bz1 );
      auto [vex1, vey1, vez1, kappa1, we2] = ExB_drift_rel( ex1, ey1, ez1, bx1, by1, bz1 );

      //-------------------------------------------------- 
      // magnetic field unit vector b at new location

      //auto [bnx1, bny1, bnz1] = mag_unit_vec(bx1, by1, bz1);
      auto [bnx1, bny1, bnz1] = mag_unit_vec_rel( ex1, ey1, ez1, bx1, by1, bz1, we2);

      //-------------------------------------------------- 
      // location update
      // NOTE: assuming mu -> 0 zero during the time step

      // Construct new gyrovelocity from conservation of magnetic moment
      // magnetic moment = m u_g^2/2 B_0 \gamma
        
      mu = 0.0; // NOTE: synchrotron losses are assumed to bring mag. mom. to zero
      real_long b1 = sqrt( bx1*bx1 + by1*by1 + bz1*bz1 );
      ug2n = mu*2.0*b1*kappa1/container.m;

      const real_long k1 = sqrt(1.0 + upar01*upar01 + ug2n);     // gamma

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
      real_long Hx = R1x - (R0x + c*vn1x);
      real_long Hy = R1y - (R0y + c*vn1y);
      real_long Hz = R1z - (R0z + c*vn1z);
      real_long H = sqrt(Hx*Hx + Hy*Hy + Hz*Hz);

      //-------------------------------------------------- 
      // guiding center location update
      if(D>=1) R1x = R0x + vn1x*c;
      if(D>=2) R1y = R0y + vn1y*c;
      if(D>=3) R1z = R0z + vn1z*c;


      //-------------------------------------------------- 
      if(false) {
      //if(n == 0) {
      //if(mu > 1.0) {
        //crash_flag = true;
          
        real_long b1 = sqrt( bx1*bx1 + by1*by1 + bz1*bz1 );
        real_long e1 = sqrt( ex1*ex1 + ey1*ey1 + ez1*ez1 );

        // E.B violation; i.e .how bad is the pure ExB drift assumption
        real_long EB0 = (ex0*bx0 + ey0*by0 + ez0*bz0)/b0/b0;
        real_long EB1 = (ex1*bx1 + ey1*by1 + ez1*bz1)/b1/b1;
        real_long bn0 = bnx0*bnx0 + bny0*bny0 + bnz0*bnz0;
        real_long bn1 = bnx1*bnx1 + bny1*bny1 + bnz1*bnz1;

        real_long egtb0 = e0/b0;
        real_long egtb1 = e1/b1;

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
        //<< " ex0:" << ex0 << " ey0:" << ey0 << " ez0:" << ez0
        //<< " bx0:" << bx0 << " by0:" << by0 << " bz0:" << bz0
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
        << "\n"
        << " mu:" << mu 
        << " ug:" << sqrt(ug2) << " ugn:" << sqrt(ug2n)
        << "\n"
        << " vex0:" << vex0 << " vey:" << vey0 << " vez:" << vez0 << " kappa0:" << kappa0
        << "\n"
        << " vex1:" << vex1 << " vey:" << vey1 << " vez:" << vez1 << " kappa1:" << kappa1
        << "\n\n";
      std::cout << std::flush;

      }

      // exit if converged
      if(H < 1e-5) break;

    }//end of iteration

    if(crash_flag) assert(false);


    vel[0][n] = static_cast<real_prtcl>( un1x );
    vel[1][n] = static_cast<real_prtcl>( un1y );
    vel[2][n] = static_cast<real_prtcl>( un1z );

    // position update from iteration, new location is following gyro center position
    if(D>=1) loc[0][n] = R1x;
    if(D>=2) loc[1][n] = R1y;
    if(D>=3) loc[2][n] = R1z;  

    // store also the field values at the new point 
    exP[n] = static_cast<real_prtcl>( ex1 );
    eyP[n] = static_cast<real_prtcl>( ey1 );
    ezP[n] = static_cast<real_prtcl>( ez1 );
    bxP[n] = static_cast<real_prtcl>( bx1 );
    byP[n] = static_cast<real_prtcl>( by1 );
    bzP[n] = static_cast<real_prtcl>( bz1 );

    bool debug_flag = 
    std::isnan(vel[0][n]) ||
    std::isnan(vel[1][n]) ||
    std::isnan(vel[2][n]) ||
    std::isnan(loc[0][n]) ||
    std::isnan(loc[1][n]) ||
    std::isnan(loc[2][n]);   

    //if(1./kinv01 > 30.0) debug_flag = true;
    if(debug_flag){
      std::cout 
        << " loc0n:" << loc[0][n] << " loc1n:" << loc[1][n] << " loc2n:" << loc[2][n]
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
  }

}



//--------------------------------------------------
// explicit template instantiation

template class pic::rGCAPusher<1,3>; // 1D3V
template class pic::rGCAPusher<2,3>; // 2D3V
template class pic::rGCAPusher<3,3>; // 3D3V

