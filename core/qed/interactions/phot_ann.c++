#include <iostream>
#include <cmath>
#include <cassert>

#include "core/qed/interactions/phot_ann.h"
#include "tools/vector.h"

#pragma GCC diagnostic push
#pragma GCC diagnostic warning "-Wunused-parameter"

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


tuple<float, float> PhotAnn::get_minmax_ene( string /*t1*/, string /*t2*/, double ene)
{
  if(ene > 0.0){

    // require e1 * e2 >= 1 and use it to get minimum required photon energy e2
    float e2 = 1.0f/ene;
    return {e2, INF};

    // NOTE: this is wrong; maximum interaction ene is always INF
    //e1 = ene;
    //e2 = 1.0/e1; 

    //if( e1 < e2 ){
    //  return {e2, e1};
    //} else {
    //  return {e1, e2};
    //}
  }

  return {0.0f, INF}; // NOTE: error branch; should not get here
}


// Exact formula for Breit-Wheeler cross section in units of sigma_T
// NOTE: sigma_T = 8 pi r2/3; this is in units of \pi r_e^2
PhotAnn::pair_float PhotAnn::comp_cross_section(
    string /*t1*/, float ux1, float uy1, float uz1,
    string /*t2*/, float ux2, float uy2, float uz2)
{

  Vec3 x1v(ux1, uy1, uz1);
  Vec3 x2v(ux2, uy2, uz2);

  float x1 = norm(x1v); // primary photon energy
  float x2 = norm(x2v); // secondary photon energy

  Vec3 om1 = x1v/x1;      // primary photon direction vector
  Vec3 om2 = x2v/x2;      // secondary photon direction vector

  float mu = dot(om1, om2); // cosine of angle between incident and target photon;
    // max(-1.0, min(dot(om1, om2), 1.0)) # cosine of angle between incident and target photons
  float xcm = sqrt(0.5*x1*x2*(1.0 - mu)); //# photon energy in COM; x_cm; g_cm
  //#bcm = sqrt(1 - 2/(x1*x2*(1-mu))) # photon COM frame; Coppi & Blandford Eq 4.3


  //if 0.5*x1*x2*(1 - mu) < 0:
  //    print('ERROR BWcs', x1, x2, om1, om2, mu)
  //    return 0.0, 1.0

  // require x_cm = g_cm > 1
  if(xcm < 1.0) return {0.0, 1.0};

  //ecm = 1.0*x;                    // e_cm = \gamma_cm
  //bcm = sqrt(1.0 - 1.0/ecm**2);   // \beta_cm

  //--------------------------------------------------
  // parameter: x = sqrt( 1 - (m_e c^2)^2/omega*omegap )
  //             = sqrt( 1 - (m_e c^2)^2/(hnu*hnup*(1-cos(mu)))
  //             = \beta_cm

  // cross section from N99; eq 64 and 59
  float x = sqrt(1.0 - 1.0/(xcm*xcm)); // \beta_cm
  float s0 = x*(1.0 - x*x)*( ( (3.0 - x*x*x*x)/(2.0*x) )*log((1.0 + x)/(1.0 - x)) - (2.0 - x*x)); 

  //s0ann = (0.25*(1/x**2)*(1-x**2))*( (3-x**4)*log((1+x)/(1-x)) + 2*x*(x**2 - 2)) #
  //s01 = 2*x*x*s0ann # eq 64 in N99; \sigma_bth = 2\beta_cm^2 \sigma_ann

  //# expression with a substitution as done by A&B book 
  //#s0 = (1/2)     *(1-x**2)* ( (3-x**4)*log( (1+x)/(1-x) ) + 2*x*(x**2-2))

  s0 *= 3.0/8.0; //# normalization factor
  //--------------------------------------------------

  //# pair creation probability dP/dtau = s_gg (1-cosa)
  float vrel = 1.0 - mu;

  return {s0, vrel};
}



void PhotAnn::interact(
  string& t1, float& ux1, float& uy1, float& uz1,
  string& t2, float& ux2, float& uy2, float& uz2) 
{

  Vec3 x1v(ux1, uy1, uz1); // four-vector of photon1
  Vec3 x2v(ux2, uy2, uz2); // four-vector of photon2

  float x1 = norm(x1v); // primary photon energy
  float x2 = norm(x2v); // secondary photon energy

  Vec3 om1 = x1v/x1;      // primary photon direction vector
  Vec3 om2 = x2v/x2;      // secondary photon direction vector

  float mu = dot(om1, om2); // cosine of angle between incident and target photon;


  // CoM frame variables; x_c
  float s  = sqrt(x1*x1 + x2*x2 + 2.0*x1*x2*mu);
  float s0 = x1 + x2;  // #sqrt(2q) in CoM
  float q  = x1*x2*(1.0 - mu);
  Vec3 svec = x1v + x2v; // or svec = x1*om1 + x2*om2

  //# CoM frame variables; x_c
  Vec3 bc = svec/(-1.0f*s0);   // lorentz transform along CoM velocity vector
  float gc = s0/sqrt(2.0*q); // lorentz factor of the boost; CoM vel
  Vec3 uv = gc*bc; // com four vel
  float v2 = dot(bc,bc); // v_c^2
  //float vc = sqrt(v2);

  //      #print(bc, gc, uv, v2, vc)

  //# boosted variables in CoM frame; x_cm 
  float gcm = sqrt(q/2.0); // prtcl energies are equal in this frame; photons also have this energy
  float vcm = sqrt(1.0 - 1.0/(gcm*gcm)); // v_cm

  // #vcm^2 = 1-2/q = 1 - 2/x1*x2*(1-mu)
  // #vcm -> 0

  //      #2 = x1*x2*(1-mu) 
  //      #2/x1*x2 = 1-mu
  //      #1 - 2/(x1*x2) = mu

  //      #1 - 2/10 = mu
  //      #mu = 1 - 0.2 = 0.8

  //      # TODO FIXME
  //      # NTOE: these are the same  angle
  //      #y = (x1 + x2*mu)/s # angle between b_cm and b_c; between photon and CoM; eq 27 in N99???
  //      #y = (s**2 - (x1 + x2)*(x2-x1))/(2*s*x1) #eq 14a in Aharonyan+ 83 for psi (substituting |k| = s

  //  build new coordinate system along svec
  Vec3<float> kvec = bc/norm(bc);            // z axis along b_c vector
  Vec3<float> jvec = unit_cross(kvec, om1);  // y axis to photon1 direction  
  Vec3<float> ivec = unit_cross(jvec, kvec); // x axis just orthogonal to others 
  Mat3<float> M(ivec, jvec, kvec);  // 3x3 rotation matrix

  float snorm = sqrt( s0*s0 - 2.0*q );
  float cosz = (x2 - x1)/snorm;  // # angle between b_c and photon

  float ymin = std::max(-1.0, (2.0 - s0)/(vcm*snorm) ); // minimum cos angle that the lepton and b_c can make
  float ymax = std::min(+1.0, (s0 - 2.0)/(vcm*snorm) ); // maximum cos angle that the lepton and b_c can make

  //--------------------------------------------------
  //  draw angles
  int niter = 0;
  float phi, xang, y, z1, z2, F;
  while(true) {
    phi = 2.0*PI*rand(); // azimuthal symmetry angle
    y = rand_ab(ymin, ymax); // # angle between b_c and lepton; lepton and CoM
                           
    // TODO mumax limit for x angle or z angle?
    // z = -1 + 2*random.rand() # angle between b_c and photon
    // y = -1 + 2*random.rand() # angle between 
    // #y = randU(-1, ymax) # angle between b_c and lepton
    //y = -1.0 + 2.0*rand(); // # angle between b_c and lepton
    xang = y*cosz + sqrt(1.0 - y*y)*sqrt(1.0 - cosz*cosz)*cos(phi); // scattering angle between photon and lepton; eq 15 in aharonyan 83

    //# DONE + or - for x experssion? DONE - according to aharonyan; + according Bottcher
    // four product scalar between electron/positron and primary/secondary photon
    z1 = (gcm*gcm)*(1.0 - vcm*xang);
    z2 = (gcm*gcm)*(1.0 + vcm*xang);

    // differential cross section angle part; F function 
    F = 0.5*( (z1/z2) + (z2/z1) + 2.0*( (1.0/z1) + (1.0/z2) ) - pow( (1.0/z1) + (1.0/z2) , 2)  );
    F *= 1.0/((1.0 + vcm)*gcm*gcm); //# normalize to [0,1]

    //# beta = vcm in N99

    //#print(z1,z2,F, x,y,z,phi)
    if( F > rand() ) break;   // accept angles
    if( niter > n_max ) break; // too many iterations
    niter += 1;
  }
  #ifdef DEBUG  
  if(niter > n_max) std::cerr << "PHOT-ANN WARNING: too many iterations" << std::endl;
  #endif

  // new primary lepton vectors in CoM frame
  //      #sinz = sqrt(1-z**2) # sin z
  //      #omr1R = array([ sinz*cos(phi), sinz*sin(phi), z ])

  float siny = sqrt(1.0 - y*y); //# sin z
  Vec3<float> omr1R( siny*cos(phi), siny*sin(phi), y ); // photon univ vector in CM frame

  //# rotate back to lab frame angles
  Mat3<float> Minv = inv( M );
  Vec3<float> omr1  = dot( Minv, omr1R ); 
  Vec3<float> omr2 = -1.0f*omr1; // other lepton has same but opposite direction (in CM frame) 

  // boost matrix back to lab frame; constructed from b_c vector
  Vec4<float> B1( gc,       -uv(0),                    -uv(1),                    -uv(2)                 );
  Vec4<float> B2(-uv(0), 1.+(gc-1.)*bc(0)*bc(0)/v2,    (gc-1.)*bc(1)*bc(0)/v2,    (gc-1.)*bc(2)*bc(0)/v2 );
  Vec4<float> B3(-uv(1),    (gc-1.)*bc(0)*bc(1)/v2, 1.+(gc-1.)*bc(1)*bc(1)/v2,    (gc-1.)*bc(2)*bc(1)/v2 );
  Vec4<float> B4(-uv(2),    (gc-1.)*bc(0)*bc(2)/v2,    (gc-1.)*bc(1)*bc(2)/v2, 1.+(gc-1.)*bc(2)*bc(2)/v2 );
  Mat4<float> B(B1, B2, B3, B4);

  // four momenta of leptons
  Vec4<float> zp1(gcm, gcm*vcm*omr1(0), gcm*vcm*omr1(1), gcm*vcm*omr1(2));
  Vec4<float> zp2(gcm, gcm*vcm*omr2(0), gcm*vcm*omr2(1), gcm*vcm*omr2(2));

  //  boost 
  Vec4<float> zpp1 = dot(B, zp1);
  Vec4<float> zpp2 = dot(B, zp2);

  float gp1 = zpp1(0); // energy of primary lepton
  float gp2 = zpp2(0); // energy of secondary lepton

  Vec3<float> ve1( zpp1(1)/gp1, zpp1(2)/gp1, zpp1(3)/gp1 );
  Vec3<float> ve2( zpp2(1)/gp2, zpp2(2)/gp2, zpp2(3)/gp2 );

  //--------------------------------------------------
  // alternatively; can deduce secondary particle from energy and momentum conservation and only use one boost
  //      #gp2 = s0 - gp1
  //      #zpp2 = svec - zpp1




  #ifdef DEBUG
  // --------------------------------------------------
  // test energy conservation
  //--------------------------------------------------
  // # test energy conservation # NOTE: we can remove these debug tests if needed
  if(true){
    float enec = x1 + x2 - (gp1 + gp2);
    Vec3<float> momc; // momc = (x1*om1 + x2*om2) - (gp1*ve1 + gp2*ve2)
    for(size_t i=0; i<3; i++) momc(i) = ( x1*om1(i) + x2*om2(i) ) - (gp1*ve1(i) + gp2*ve2(i) );
    float moms = sum(momc);

    bool ts[8]; // t1,t2,t3,t4,t5,t6,t7,t8 = False,False,False,False,False,False,False,False
    for(size_t i = 0; i < 2; i++) ts[i] = false;

    float nom1  = norm(om1);
    float nom2  = norm(om2);

    if(abs(enec)         > tol) ts[0] = true;
    if(abs(sum(momc)) > tol) ts[1] = true;

    if(ts[0] ||
       ts[1] ) { 

      std::cerr << "ERROR PHOT-ANN:" << std::endl;
      for(size_t i = 0; i < 2; i++) { std::cerr << i << " " << ts[i] << std::endl; }

      std::cerr << "x1v, x2v  " <<  x1v   << " " <<  x2v << std::endl;
      std::cerr << "gm, gp    " <<  gp1   << " " <<  gp2  << std::endl;
      std::cerr << "mu,s0,s,q " <<  mu    << " " <<  s0    << " " <<  s << " " << q  << std::endl;
      std::cerr << "x,x1,enec " <<  x1    << " " <<  x2    << " " <<  enec << std::endl;
      std::cerr << "momc      " <<  moms  << " " <<  momc  << std::endl;
      std::cerr << "|om|,|om1|" <<  nom1  << " " <<  nom2  << std::endl;
      std::cerr << "|ve| |vp| " <<  ve1   << " " <<  ve2 << std::endl;
      std::cerr << "th,phi,om " <<  siny  << " " <<  phi   << "om1 " << om1 << "om2 " <<  om2  << std::endl;
      std::cerr << "zp1,zp2   " <<  zp1   << " " <<  zp2   << std::endl;
      std::cerr << "B         " <<  B      <<                 std::endl;
      std::cerr << "omr1,omr2 " <<  omr1  << " " <<  omr2  << std::endl;
    }
  }
  #endif

  // is this flip needed? seems so; makes routine independent of e-/e+
  if(rand() < 0.5) {
    t1 = "e-";
    t2 = "e+";
  } else {
    t1 = "e+";
    t2 = "e-";
  }

  ux1 = zpp1(1);
  uy1 = zpp1(2);
  uz1 = zpp1(3);

  ux2 = zpp2(1);
  uy2 = zpp2(2);
  uz2 = zpp2(3);

  return;
}

} // end of namespace qed

#pragma GCC diagnostic pop
