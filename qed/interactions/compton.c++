#include <iostream>
#include <cmath>
#include <cassert>
#include "compton.h"
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


tuple<float_p, float_p> Compton::get_minmax_ene( string t1, string t2, double ene)
{
  return {0.0f, INF};


  // TODO brainstorming here below; routine does not seem possible
  // skip everything under minimum energy (depending on type)
  //if( (t1 == "e-" || t1 == "e+") && ene < minz ) return {INF, INF};
  //if( (t1 == "ph"              ) && ene < minx ) return {INF, INF};
  //// check also target and make it mininum limit
  //if(t2 == "e-" || t2 == "e+") return {minz, INF}; // pair
  //if(t2 == "ph")               return {minx, INF}; // photon

  //// never get here
  //assert(false); 
}

Compton::pair_float Compton::comp_cross_section(
    string t1, float_p ux1, float_p uy1, float_p uz1,
    string t2, float_p ux2, float_p uy2, float_p uz2)
{

  Vec3<float_p> zv, xv;
  if( (t1 == "e-" || t1 == "e+") && (t2 == "ph") ) {
    zv.set(ux1, uy1, uz1); 
    xv.set(ux2, uy2, uz2); 
  } else if( (t2 == "e-" || t2 == "e+") && (t1 == "ph") ) {
    zv.set(ux2, uy2, uz2); 
    xv.set(ux1, uy1, uz1); 
  } else {
    // incompatible particle types
    assert(false);
  }

  float_p x = norm(xv);                   // photon energy
  Vec3<float_p> om = xv/x;                // photon direction vector \omega
  float_p gam = sqrt(1.0 + dot(zv,zv));   // prtcl gamma
  Vec3<float_p> beta = zv/gam;            // prtcl 3-velocity \beta
  float_p beta0 = norm(beta);
  float_p mu = dot(om, beta)/beta0;       // cosine of angle between incident photon and electron
  float_p s = 0.5*x*gam*(1.0 - beta0*mu); // relativistic invariant


  //# Approximate Taylor expansion for Compton total 
  // cross-section if xi<0.01, error about 5e-6
  float_p s0 = 0.0; 
  if(s < 0.01) {
    s0 = 1.0 - 2.0*s + 5.2*s*s - 9.1*s*s*s;

    //# higher-order terms, not needed if xi<0.01 
    //# + 1144.0 * xi**4 / 35.0  - 544.0 * xi**5 / 7. 
    //#+ 1892.0 * xi**6 / 21.0   
  } else {
    // Exact formula for Klein-Nishina cross section in units of sigma_T
    s0 = (1.0/(s*s))*(4.0 + (s - 2.0 - 2.0/s)*log(1.0 + 2.0*s) + 2.0*s*s*(1.0 + s)/pow(1.0 + 2.0*s, 2) );
    s0 *= 3.0/8.0;  //change from \sigma_0 to \sigma_T
  }

  // compton scattering of two particles will have relative velocity:
  //# vrel = 1-zeta = 1 - \beta mu so that dP/dtau ~ s_0 *(1-beta mu)/2
  //#vrel =  (1 - beta0*mu)/gam/x # TODO check if correct: normalize number densities with lorentz factors
  float_p fkin = (1.0 - beta0*mu);

  //std::cout<< "comp s:" << s << " x:" << x << " s0:" << s0 << " beta0:" << beta0 << " mu:" << mu << std::endl;

  return {s0, fkin};
}


tuple<float_p, float_p> Compton::accumulate(
    string t1, float_p e1, string t2, float_p e2)
{
  if( (t1 == "e-" || t1 == "e+") && e1 > ming) return {1.0f, 1.0f}; // do not accumulate rel prtcl
  if( (t2 == "e-" || t2 == "e+") && e2 > ming) return {1.0f, 1.0f}; // do not accumulate rel prtcl


  // accumulation factor; forces electron energy changes to be ~0.1
  float_p f = std::max(1.0f, 0.01f/(e1*e2));

  f = std::min(1.0e3f, f); // cap f 

  float_p f1 = t1 == "ph" ? f : 1.0f;
  float_p f2 = t2 == "ph" ? f : 1.0f;

  return {f1,f2};
}
  
void Compton::interact(
  string& t1, float_p& ux1, float_p& uy1, float_p& uz1,
  string& t2, float_p& ux2, float_p& uy2, float_p& uz2) 
{

  //--------------------------------------------------
  Vec3<float_p> zv, xv;
  if( (t1 == "e-" || t1 == "e+") && (t2 == "ph") ) {
    zv.set(ux1, uy1, uz1); 
    xv.set(ux2, uy2, uz2); 
  } else if( (t2 == "e-" || t2 == "e+") && (t1 == "ph") ) {
    zv.set(ux2, uy2, uz2); 
    xv.set(ux1, uy1, uz1); 
  } else {
    // incompatible particle types
    assert(false);
  }

  //--------------------------------------------------

  float_p x0 = norm(xv);                  // photon energy
  Vec3<float_p> om0 = xv/x0;              // photon direction vector \omega
  float_p gam0 = sqrt(1.0f + dot(zv,zv));  // prtcl gamma
  Vec3<float_p> beta0 = zv/gam0;          // prtcl 3-velocity \beta
  Vec3<float_p> bdir = beta0/norm(beta0); // electron direction vector


  //#--------------------------------------------------
  // check for Thomson regime
  //float_p gam_min = sqrt(1.0f + minz*minz); // min limit from class member
  bool flag_thomson = false;
  if( gam0 < ming && gam0*x0 < minx2z ) flag_thomson = true;


  //#--------------------------------------------------
  // boost to electron mom frame

  // NOTE: renaming variables so that B matrix matches with other processes
  float_p gc = gam0;               //# frame lorentz factor
  Vec3<float_p> uv = gam0*beta0;  //# momentum 3-vec of the farme
  Vec3<float_p> bc = beta0;             //# 3-vel vec of the frame
  float_p v2 = dot(bc,bc);         //# v_c^2
  //Vec3<float_p> beta_dir = beta0/norm(beta0); // vector into direction of electron velocity FIXME
          
  //# boost matrix 
  Vec4<float_p> B1( gc,       -uv(0),                    -uv(1),                    -uv(2)                 );
  Vec4<float_p> B2(-uv(0), 1.+(gc-1.)*bc(0)*bc(0)/v2,    (gc-1.)*bc(1)*bc(0)/v2,    (gc-1.)*bc(2)*bc(0)/v2 );
  Vec4<float_p> B3(-uv(1),    (gc-1.)*bc(0)*bc(1)/v2, 1.+(gc-1.)*bc(1)*bc(1)/v2,    (gc-1.)*bc(2)*bc(1)/v2 );
  Vec4<float_p> B4(-uv(2),    (gc-1.)*bc(0)*bc(2)/v2,    (gc-1.)*bc(1)*bc(2)/v2, 1.+(gc-1.)*bc(2)*bc(2)/v2 );
  Mat4<float_p> B(B1, B2, B3, B4);

  //--------------------------------------------------

  Vec4<float_p> k0( x0, x0*om0(0), x0*om0(1), x0*om0(2) );  //# 4-mom of the photon
  Vec4<float_p> k0_R = dot(B, k0);  //# boosted 4-mom of photon in electron rest frame

  float_p x0_R = k0_R(0); //# energy of photon in R frame
  Vec3<float_p> om_R(k0_R(1)/x0_R, k0_R(2)/x0_R, k0_R(3)/x0_R);  //# direction of photon in R frame

  // --------------------------------------------------
  // test that electron becomes at rest after boost
  // z0 = gam0*np.array([1, beta0[0], beta0[1], beta0[2]])  # 4-mom of the photon
  // z0_R = np.dot(B, z0) # boosted 4-mom of photon in electron rest frame
  //  DONE; confirmed to work

  //--------------------------------------------------
  // choose scattering coordinates (i,j,k); along photon direction
  Vec3<float_p> kvec = om_R;                       // along photon direction
  //Vec3<float_p> jvec = unit_cross(kvec, beta0); // k x v; incoming electron is in the (i,k) plane FIXME
  Vec3<float_p> jvec = unit_cross(kvec, bdir); // k x v; incoming electron is in the (i,k) plane FIXME
  Vec3<float_p> ivec = unit_cross(jvec, kvec);     // # complement the right-handed basis

  // transformation matrix between lab frame and ijk frame
  Mat3<float_p> M(ivec, jvec, kvec);           // 3x3 rotation matrix

  //--------------------------------------------------
  //scatter until physically realistic event is found
  int niter = 0;
  float_p phi_R, mu_R, x1_R, F;
  while(true) {
    phi_R = 2.0f*PI*rand(); // candidate symmetry/azimuth angle in R frame
    mu_R = -1.0f + 2.0f*rand(); // candidate latitude/scattering angle between incident and outg photon

    x1_R = x0_R/(1.0f + x0_R*(1.0f - mu_R)); // scattered photon energy in R frame

    // angle dependent part of differential cross section
    F = 0.5f*pow(x1_R/x0_R, 2)*(-1.0f + (x1_R/x0_R) + (x0_R/x1_R) + mu_R*mu_R );

    if( F > rand() ) break;   // accept angles
    if( niter > 10000 ) break; // too many iterations
    niter += 1;
  }
        
  if(niter > 10000) std::cout << "COMPTON WARNING: too many iterations" << std::endl;


  //# construct new photon vector based on the angles
  //      #om1_Rijk  = np.array([mu, sinth*np.cos(phi), sinth*np.sin(phi)])
  float_p sinz = sqrt(1.0f - mu_R*mu_R);
  Vec3<float_p> om1_Rijk( sinz*sin(phi_R), sinz*cos(phi_R), mu_R);

  //# rotate back to original axis
  Mat3<float_p> Minv = inv( M );
  Vec3<float_p> om1_R = dot( Minv, om1_Rijk );

  //#-------------------------------------------------- 
  //# transform back to lab frame


  // boost matrix back to lab frame; constructed from -b_c vector
  //float_p gc = gam0; // # frame lorentz factor
  uv = -1.0f*uv;  // flip momentum 3-vec of the farme
  bc = -1.0f*bc;  //  3-vel vec of the frame
  //v2 = np.dot(bc,bc) # v_c^2
    
  //# boost matrix 
  Vec4<float_p> G1( gc,       -uv(0),                    -uv(1),                    -uv(2)                 );
  Vec4<float_p> G2(-uv(0), 1.+(gc-1.)*bc(0)*bc(0)/v2,    (gc-1.)*bc(1)*bc(0)/v2,    (gc-1.)*bc(2)*bc(0)/v2 );
  Vec4<float_p> G3(-uv(1),    (gc-1.)*bc(0)*bc(1)/v2, 1.+(gc-1.)*bc(1)*bc(1)/v2,    (gc-1.)*bc(2)*bc(1)/v2 );
  Vec4<float_p> G4(-uv(2),    (gc-1.)*bc(0)*bc(2)/v2,    (gc-1.)*bc(1)*bc(2)/v2, 1.+(gc-1.)*bc(2)*bc(2)/v2 );
  Mat4<float_p> G(G1, G2, G3, G4);

  Vec4<float_p> k1_R( x1_R, x1_R*om1_R(0), x1_R*om1_R(1), x1_R*om1_R(2) ); //# construct photon 4-mom to be boosted
  Vec4<float_p> k1 = dot(G, k1_R); // # de-boost back to lab frame


  //# --------------------------------------------------
  //# construct new scattered quantities
  float_p x1 = k1(0);
  Vec3<float_p>  om1( k1(1)/x1, k1(2)/x1, k1(3)/x1 );

  // accumulation factor; NOTE: t1/t2 order does not matter here
  auto [facc1, facc2] = accumulate(t1, gam0, t2, x0);
  facc1 = do_accumulate ? facc1 : 1.0f;
  facc2 = do_accumulate ? facc2 : 1.0f;
  float_p facc_in = t1 == "ph" ? facc1 : facc2; // pick photon as the accumulated quantity
  //facc = 1.0f; // FIXME never accumulate losses

  // every energy transcation is is enhanced by this factor
  //facc *= wtar2wini;

  // limit change to 0.1\gamma_0
  //facc = std::min(facc, (0.3f*gam0 - x0)/x1);
  //x1 *= facc;

  // limit facc so that gam1 is > 1
  // gam1 = gam0 + x0 - facc*x1; 
  // 1 > gam0 + x0 - facc*x1
  // 1 -gam0 - x0 > - facc*x1
  // gam0 + x0 - 1> facc*x1
  // (gam0 + x0 - 1)/x1 > facc
  float_p facc_max = (gam0 + x0 - (1.0f + 1e-4f) )/x1;
  float_p facc = std::min(facc_in, facc_max );

  //minimum limit
  facc = std::max(1.0f, facc); 

  // TODO which one is right?
  //# scattered electon variables from ene and mom conservation
  //float_p gam1 = gam0 + facc*(x0 - x1); // ver1
  float_p gam1 = gam0 + x0 - facc*x1; // ver2
  //float_p gam1 = gam0 + x0 - x1; // ver3// facc in x1

  //gam1 = gam0 + facc*(x0 - x1) 
  //     = gam0 + facc*x0 - facc*x1 
  //     = gam0 + 0.1*x0/x0*gam0 - 0.1*x1/x0*gam0
  //     = gam0 + 0.1*x0/x0*gam0 - 0.1*x1/x0*gam0
  //     = gam0 + 0.1/gam0 - 0.1*x1/x0*gam0
  //     = gam0 +~0.1 -~ 0.1*gam0/x0*gam0
  //     = gam0 +~0.1 -~ 0.1/x0


  //Vec3<float_p> beta1 = (gam0*beta0 + x0*om0 - x1*om1)/gam1;
  Vec3<float_p> beta1(0.0f, 0.0f, 0.0f);
  for(size_t i=0; i<3; i++) beta1(i) = (gam0*beta0(i) + (x0*om0(i) - facc*x1*om1(i)) )/gam1; // FIXME
  //for(size_t i=0; i<3; i++) beta1(i) = (gam0*beta0(i) + (x0*om0(i) - x1*om1(i)) )/gam1;


  if( (t1 == "e-" || t1 == "e+") && (t2 == "ph") ) {

    if(! no_electron_update) {
      ux1 = gam1*beta1(0);
      uy1 = gam1*beta1(1);
      uz1 = gam1*beta1(2);
    }

    if(! no_photon_update) {
      ux2 = x1*om1(0);
      uy2 = x1*om1(1);
      uz2 = x1*om1(2);
    }

  } else if( (t2 == "e-" || t2 == "e+") && (t1 == "ph") ) {

    if(! no_photon_update) {
      ux1 = x1*om1(0);
      uy1 = x1*om1(1);
      uz1 = x1*om1(2);
    }

    if(! no_electron_update) {
      ux2 = gam1*beta1(0);
      uy2 = gam1*beta1(1);
      uz2 = gam1*beta1(2);
    }
  }


  // --------------------------------------------------
  // test energy conservation
  //--------------------------------------------------
  // # test energy conservation # NOTE: we can remove these debug tests if needed
  //if(true && !(do_accumulate) ){
  //if(false){
  if(true){

    float_p enec = gam1 + x1 - (x0 + gam0);

    Vec3<float_p> momc; 
    for(size_t i=0; i<3; i++) momc(i) = x1*om1(i) + gam1*beta1(i)  - (x0*om0(i) + gam0*beta0(i) );
    float_p moms = sum(momc);

    bool ts[8]; // t1,t2,t3,t4,t5,t6,t7,t8 = False,False,False,False,False,False,False,False
    for(size_t i = 0; i < 5; i++) ts[i] = false;

    float_p tol = 3.0e-5;

    float_p nom0  = norm(om0);
    float_p nom1  = norm(om1);

    if(abs(enec)      > tol) ts[0] = true;
    if(abs(sum(momc)) > tol) ts[1] = true;
    if(gam0 < 1.0)           ts[2] = true;
    if(gam1 < 1.0)           ts[3] = true;
    //if(flag_thomson)         ts[4] = true; // debug Thomson regime


    if(ts[0] ||
       ts[1] ||
       ts[2] ||
       ts[3] ||
       ts[4] ) { 

      std::cout << "ERROR COMPTON:" << std::endl;
      for(size_t i = 0; i < 5; i++) { std::cout << i << " " << ts[i] << std::endl; }

      std::cout << "x0v, x1v  " <<  xv    << " " <<  k1   << std::endl;
      std::cout << "g0, g1    " <<  gam0  << " " <<  gam1  << std::endl;
      //std::cout << "mu,s0,s,q " <<  mu_R  << " " <<  s0    << " " <<  s << " " << q  << std::endl;
      std::cout << "x,x1,enec " <<  x0    << " " <<  x1    << " " <<  enec << std::endl;
      std::cout << "momc      " <<  moms  << " " <<  momc  << std::endl;
      std::cout << "|om0||om1|" <<  nom0  << " " <<  nom1  << std::endl;
      std::cout << "facc      " <<  facc  << " de" <<  x0-x1  << std::endl;
      std::cout << "facc lim  " <<  (0.5f*gam0 -x0)/x1  << std::endl;
      std::cout << "facc max  " <<  facc_max  << std::endl;
      std::cout << "facc in  "  <<  facc_in  << std::endl;
    }
  }


  return;
}


} // end of namespace qed
