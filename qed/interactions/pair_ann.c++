#include <iostream>
#include <cmath>
#include "pair_ann.h"
#include "../../tools/vector.h"


namespace qed {
  using std::string;
  using std::tuple;
  using std::sqrt;

  using toolbox::Vec3;
  using toolbox::Mat3;
  using toolbox::Vec4;
  using toolbox::Mat4;
  using toolbox::norm;
  using toolbox::dot;
  using toolbox::cross;
  using toolbox::unit_cross;
  using toolbox::sum;
  using toolbox::inv;


tuple<double, double> PairAnn::get_minmax_ene( string t1, string t2)
{
  return {0.0, 3.0};
}


double PairAnn::comp_cross_section(
    string t1, double ux1, double uy1, double uz1,
    string t2, double ux2, double uy2, double uz2)
{

  double zp = norm(ux1, uy1, uz1); // z_+
  double zm = norm(ux2, uy2, uz2); // z_-
  double gamp = sqrt(zp*zp + 1.0); // \gamma_+
  double gamm = sqrt(zm*zm + 1.0); // \gamma_-

  double zeta = ( ux1*ux2 + uy1*uy2 + uz1*uz2 )/( zm*zp ); //#angle between pair's momentum
  
  //qe = max(1+EPS, gamp*gamm - zp*zm*zeta) 
  double qe   =  gamp*gamm - zp*zm*zeta; //# product of 4-moms; q_e = z_- . z_+; = 2gamma_cm^2 - 1 = gamma_r
  double q = qe + 1.0;          //# q_e = q + 1
  double bcm = sqrt( (q-2.0)/q ); //# \beta_cm; matches beta' in Coppi & Blandford
  //double gcm = 1.0/sqrt(1.0 - bcm*bcm); //# gamma_cm

  //# expression via relativistic invariant x = sqrt( p_- . p_+ )
  double x = bcm;                  //# free variable
  double s0 = (0.25*(1.0/x*x)*(1.0-x*x))*( (3.0-x*x*x*x)*log((1.0+x)/(1.0-x)) + 2.0*x*(x*x - 2.0)); //# ver2
  s0 *= 3.0/8.0; //# normalization to rates; mathches with Coppi & Blandford eq 3.1

  //# ver2; matches above
  //#s0 = 1/(4*bcm*gcm**2)*( (1/bcm)*(2 + 2/gcm**2 - 1/gcm**4)*np.log( (1+bcm)/(1-bcm) ) - 2 - 2/gcm**2)
  //#s0 *= 3.0/8.0 # normalization to rates; mathches with Coppi & Blandford eq 3.1
    

  //#--------------------------------------------------
  //# relative velocity calculation

  //# testing which one is correct.... and none of these are correct.
  //# 1 - ( beta_+ \beta_- cos(\theta) )
  //#vrel = 1 - (zp/gamp)*(zm/gamm)*zeta # seems most correct according to qed3 test results
  //#vrel = np.sqrt( (zp*zm*zeta)**2 - 1/(gamp*gamm)) # produces hi-ene peak; wrong
  //#vrel = (zp/gamp)*(zm/gamm)*zeta 
  //# NOTE: all of these are equal up to numerical precision; last one is roughly fastest

  //# fkin; eq 3.3 in Coppi & Blandford
  //#betam = zm/gamm
  //#betap = zp/gamp
  //#fkin1 = np.sqrt( betam**2 + betap**2 - (betam*betap)**2*(1-zeta**2) - 2*betam*betap*zeta )

  //# kinematic factor based on Svensson 82 Eq 9
  //#gamcm = 1/np.sqrt(1 - bcm**2) # gamma_cm
  //#fkin2 = 2*bcm*gamcm**2/(gamm*gamp)
  //#fkin2 = bcm*q/(gamm*gamp) # equal expression to above

  // kinematic factor based on Svensson 82 Eq 8; qe = q_+ . q_-
  double fkin = sqrt( qe*qe - 1.0 )/(gamp*gamm);

  //if qe**2 -1 < 0: print('pair-ann cs:', s0, fkin, qe, gamp, gamm, zeta, x)
    
  return s0*fkin;
}


//tuple<
//    string, double, double, double,
//    string, double, double, double>
//PairAnn::interact(
//  string t1, double ux1, double uy1, double uz1,
//  string t2, double ux2, double uy2, double uz2) 
//{
//
//  //      particle 1          particle 2
//  return {t3, ux3, uy3, uz3,  t4, ux4, uy4, uz4};
//}
  
void PairAnn::interact(
  string& t1, double& ux1, double& uy1, double& uz1,
  string& t2, double& ux2, double& uy2, double& uz2) 
{
  Vec3<double> zmvec(ux1, uy1, uz1);
  Vec3<double> zpvec(ux2, uy2, uz2);

  double zm = norm(zmvec); //# electron momenta z_-
  double zp = norm(zpvec); //# positron momenta z_+

  double gamm = sqrt(zm*zm + 1.0); //# electron gamma_-
  double gamp = sqrt(zp*zp + 1.0); //# positron gamma_+

  Vec3<double> omm = zmvec/zm; //# direction of electron momenta Omega_-
  Vec3<double> omp = zpvec/zp; //# direction of positron momenta Omega_+
        
  double zeta = dot(omm, omp); //#angle between electron and positron momenta

  double s0 = gamm + gamp;                        //#s0 = x + x1 #sqrt(2q) in CoM
  double s  = sqrt(zm*zm + zp*zp + 2*zm*zp*zeta); //#s  = sqrt(x**2  + x1**2 + 2*x*x1*mu)
  double q  = gamp*gamm - zp*zm*zeta + 1.0;       //#q  = x*x1*(1-mu)
  //      #q = s0**2 - s**2

  //Vec3 svec = zp*omp + zm*omm;                  //#svec = x*om + x1*om1
  Vec3<double> svec_a = zp*omp;                  //#svec = x*om + x1*om1
  Vec3<double> svec_b = zm*omm;                  //#svec = x*om + x1*om1
  Vec3<double> svec = svec_a + svec_b;

  //# CoM frame variables; x_c
  Vec3<double> bc = svec/(-1.0*s0); //#lorentz transform along CoM velocity vector
  double gc = s0/sqrt(2.0*q); //# lorentz factor of the boost; CoM vel
  Vec3<double> uv = gc*bc; //# com four vel
  double v2 = dot(bc,bc); //# v_c^2
  double vc = sqrt(v2);

  //--------------------------------------------------
  //# boosted variables in CoM frame; x_cm 
  double gcm = sqrt(q/2.0); // # prtcl energies are equal in this frame; photons also have this energy
  double vcm = sqrt(1.0-1.0/(gcm*gcm)); //# v_cm

  //# angle between b_cm and b_c; electron and CoM 
  double y = (1.0/vc/vcm)*((gamp - gamm)/(gamp + gamm));

  //# build new coordinate system along svec
  Vec3<double> kvec = bc/norm(bc); //# z axis along b_c vector
  Vec3<double> jvec = unit_cross(kvec, omm); //# y axis to electron direction  # np.cross(omp, omm)/sqrt(1 - zeta**2)
  Vec3<double> ivec = unit_cross(jvec, kvec); // # x axis just orthogonal to others 
                                      // # ( (zp + zm*zeta)*omm - (zm + zp*zeta)*omp )/(s*sqrt(1 - zeta**2))
  Mat3<double> M(ivec, jvec, kvec); // 3x3 rotation vector


  //#--------------------------------------------------
  //# draw angles
  int niter = 0;

  double cosz, phi, xang, z1, z2, F; // temp variables
  while(true) {
    //# angle between k_cm and b_c; photon and CoM
    cosz= -1.0 + 2.0*rand();
    phi = 2.0*PI*rand();

    //# angle between k_cm and b_cm; photon and electron
    xang = y*cosz + sqrt(1.0 - y*y)*sqrt(1.0 - cosz*cosz)*cos(phi);

    //# four product scalar between electron/positron and primary/secondary photon
    z1 = (gcm*gcm)*(1.0 - vcm*xang);
    z2 = (gcm*gcm)*(1.0 + vcm*xang);

    //# differential cross section angle part; F function 
    F = 0.5*( (z1/z2) + (z2/z1) + 2.0*( (1.0/z1) + (1.0/z2) ) - pow( (1.0/z1) + (1.0/z2), 2)  );
    F *= 1.0/((1.0 + vcm)*gcm*gcm); //# normalize to [0,1]

    if( F > rand() ) break; // accept angles
    if( niter > 1000 ) break; // too many iterations
    niter += 1;
  }

  //# new photon vectors in CoM frame
  double sinz = sqrt(1.0 - cosz*cosz); //# sin z
  Vec3<double> omrR( sinz*cos(phi), sinz*sin(phi), cosz );

  //# rotate back to lab frame angles
  Mat3<double> Minv = inv( M );
  Vec3<double> omr0 = dot( Minv, omrR );
  Vec3<double> omr1 = -1.0*omr0; //# other photon has same but opposite direction # np.dot( np.linalg.inv(M), omr1R ) 


  //# boost matrix back to lab frame; constructed from b_c vector
  Vec4<double> B1( gc,       -uv(0),                    -uv(1),                    -uv(2)                 );
  Vec4<double> B2(-uv(0), 1.+(gc-1.)*bc(0)*bc(0)/v2,    (gc-1.)*bc(1)*bc(0)/v2,    (gc-1.)*bc(2)*bc(0)/v2 );
  Vec4<double> B3(-uv(1),    (gc-1.)*bc(0)*bc(1)/v2, 1.+(gc-1.)*bc(1)*bc(1)/v2,    (gc-1.)*bc(2)*bc(1)/v2 );
  Vec4<double> B4(-uv(2),    (gc-1.)*bc(0)*bc(2)/v2,    (gc-1.)*bc(1)*bc(2)/v2, 1.+(gc-1.)*bc(2)*bc(2)/v2 );
  Mat4<double> B(B1, B2, B3, B4);

  //# four momenta of photons
  Vec4<double> xp0(gcm, gcm*omr0(0), gcm*omr0(1), gcm*omr0(2));
  Vec4<double> xp1(gcm, gcm*omr1(0), gcm*omr1(1), gcm*omr1(2));

  //# boost 
  Vec4<double> xpp0 = dot(B, xp0);
  Vec4<double> xpp1 = dot(B, xp1);

  double x0 = xpp0(0);  //# energy of primary photon
  double x1 = xpp1(0); //# energy of secondary photon

  Vec3<double> om0( xpp0(1), xpp0(2), xpp0(3)); ///x
  Vec3<double> om1( xpp1(1), xpp1(2), xpp1(3)); ///x1

  om0 = om0/x0;
  om1 = om1/x1;

  //--------------------------------------------------
  // # test energy conservation # NOTE: we can remove these debug tests if needed
  if(true){
    double enec = gamm + gamp - (x0 + x1);
    Vec3<double> momc; // zmvec + zpvec; - (x*om + x1*om1);
    for(size_t i=0; i<3; i++) momc(i) = zmvec(i) + zpvec(i) - ( x0*om0(i) + x1*om1(i) );
    double moms = sum(momc);

    double nom1i = norm(omm);
    double nom2i = norm(omp);
    double nom1  = norm(om0);
    double nom2  = norm(om1);

    bool ts[8]; // t1,t2,t3,t4,t5,t6,t7,t8 = False,False,False,False,False,False,False,False
    for(size_t i=0; i<8; i++)    ts[i] = false; 

    if(abs(enec) > 1.0e-12)      ts[0] = true; 
    if(x0 < 0.0)                 ts[1] = true;
    if(x1 < 0.0)                 ts[2] = true;
    if(abs(nom1i)-1. > 1.0e-12)  ts[3] = true;
    if(abs(nom2i)-1. > 1.0e-12)  ts[4] = true;
    if(abs(nom1)-1. > 1.0e-12)   ts[5] = true;
    if(abs(nom2)-1. > 1.0e-12)   ts[6] = true;
    if(abs(moms) > 1.0e-12)      ts[7] = true;

    //if(true) {
    if(ts[0] ||
       ts[1] ||
       ts[2] ||
       ts[3] ||
       ts[4] ||
       ts[5] ||
       ts[6] ||
       ts[7] ) {
      std::cout << "ERROR PAIR-ANN:" << std::endl;
      for(size_t i = 0; i < 8; i++) { std::cout << i << " " << ts[i] << std::endl; }
    
      std::cout << "zm, zp    " <<  zmvec << " " <<  zpvec << std::endl;
      std::cout << "gm, gp    " <<  gamm  << " " <<  gamp  << std::endl;
      std::cout << "z,s0,s,q  " <<  zeta  << " " <<  s0    << " " <<  s << " " << q  << std::endl;
      std::cout << "x,x1,enec " <<  x0    << " " <<  x1    << " " <<  enec << std::endl;
      std::cout << "|omi||omi|" <<  nom1i << " " <<  nom2i << std::endl;
      std::cout << "|om|,|om1|" <<  nom1  << " " <<  nom2  << std::endl;
      std::cout << "momc      " <<  moms  << " " <<  momc  << std::endl;
      std::cout << "th,phi,om " <<  cosz  << " " <<  phi   << "om " << om0 << "om1 " <<  om1  << std::endl;
      std::cout << "xpp0,xpp1 " <<  xpp0  << " " <<  xpp1  << std::endl;
      std::cout << "B         " <<  B      <<                 std::endl;
      std::cout << "omr0,omr1 " <<  omr0  << " " <<  omr1  << std::endl;
      std::cout << "xp0, xp1  " <<  xp0   << " " <<  xp1   << std::endl;
    }
  }

  //# NOTE flip randomly; does not have an effect because photons are identical
  //# alternate the role of primary/secondary; otherwise algorithm is deterministic
  //#if np.random.rand() < 0.5:
  //#    x, x1  = x1,x
  //#    om,om1 = om1, om

  t1 = "ph";
  ux1 = xpp0(1);
  uy1 = xpp0(2);
  uz1 = xpp0(3);

  t2 = "ph";
  ux2 = xpp1(1);
  uy2 = xpp1(2);
  uz2 = xpp1(3);

  return;
}


} // end of namespace qed
