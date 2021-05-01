#pragma once

// Particle weight functions


// 1st order cloud-in-cell shape (CIC)
inline void W1st(double d, double* coeff){
  coeff[0] = 0.0;   // W1_im1
  coeff[1] = 1.0-d; // W1_i
  coeff[2] = d;     // W1_ip1
}


// Triangular shaped cloud;
// 2nd order quadratic spline
/*      -
*       |  3/4 - x^2                  if |x|<1/2
* W(x)=<   1/2*(3/2 - |x|)^2          if 1/2<=|x|<3/2
*       |  0                          otherwise
*       -
*/
// TODO: c0 and c2 appear wrong; should be 1/8 (3-2x)^2 ?
inline void W2nd(double d, double* coeff){

  // Sokolov higher
  // 0.5*(1-d)^2 ok
  // 3/4 - (0.5 - d)^2
  // d^2

  //NOTE: d at wings includes +-1 so W2_mp1 -> 0.5*(3/2 +- d)^2
  coeff[0] = 0.50*(0.5 - d)*(0.5 - d); //W2_im1 
  coeff[1] = 0.75 - d*d;               //W2_i   
  coeff[2] = 0.50*(0.5 + d)*(0.5 + d); //W2_ip1 
}


// 3rd oder cubic piewce spline  
inline void W3rd(double d, double* coeff){
  double d2 = d*d;
  double d3 = d*d*d;

  // OLD version
  ///coeff[0] = ( 1.0f - d3 )*1.0f/6.0f - 0.5f*( d-d2 );
  ///coeff[1] = 2.0f/3.0f - xi2 + 0.5f*d3;
  ///coeff[2] = 1.0f/6.0f + 0.5f*( d+d2-d3 );
  ///coeff[3] = d3*1.0f/6.0f;


  // Sokolov weights
  //(1-d)**3/6  == 1/6 - d/2 + d^2/2 - d^3/6        ok
  //(2-d)**3/6 - 2(1-d)**3/3 == 0.5*d^3 - d^2 + 2/3 ok
  //(1+d)**3/6 - 2d**3/3 == 1/6 - (d + d^2 - d^3)/2 ok
  //d**3/6                                          ok

  // 1/6*(4 - 6*x^2 + 3*|x|^3)   if 0<=|x|<1
  // 1/6*(2 - |x|)^3             if 1<=|x|<2
  // 0                           otherwise

  coeff[0] = 0.0;                             // W3_im2
  coeff[1] = ( 1.0 - d3 )/6.0 - 0.5*( d-d2 ); // W3_im1  (1-d^3)/6 - (d-d^3)/2
  coeff[2] = 2.0/3.0 - d2 + 0.5*d3;           // W3_i    2/3 - d^2 + d^3/2
  coeff[3] = 1.0/6.0 + 0.5*( d+d2-d3 );       // W3_ip1  1/6 - (d + d^2 - d^3)/2
  coeff[4] = d3/6.0;                          // W3_ip2  d^3/6

}

// 4th order quartic piewce spline  
// TODO: not tested
inline void W4th(double d, double* coeff){
  double d2 = d*d;
  double d3 = d*d*d;
  double d4 = d*d*d*d;

  // 115/192 + x^2 * (-5/8 + 1/4 * x^2)                          if -1/2 < x < 1/2
  // 1/96 * (55 + 4 * x * (5 - 2 * x * (15 + 2 * x * (-5 + x)))) if 1/2 <= |x| < 3/2
  // 1/384 * (5 - 2 * x)^4                                       if 3/2 <= |x| < 5/2

  coeff[0] = 1.0/384.0   - d/48.0       + d2/16.0    - d3/12.0 + d4/24.0; // W4_im2
  coeff[1] = 19.0/96.0   - d*11.0/24.0  + d2/4.0     + d3/6.0  - d4/6.0; // W4_im1
  coeff[2] = 115.0/192.0                - d2*5.0/8.0           + d4/4.0; // W4_i  //ok
  coeff[3] = 19.0/96.0   + d*11.0/24.0  + d2/4.0     - d3/6.0  - d4/6.0; // W4_ip1
  coeff[4] = 1.0/384.0   + d/48.0       + d2/16.0    + d3/12.0 + d4/24.0; // W4_ip2

}


