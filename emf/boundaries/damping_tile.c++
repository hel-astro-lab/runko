#include "damping_tile.h"


using std::exp;


// deposit current with some resistivity
//template<
//  std::size_t D, 
//  int S
//>
//void emf::damping::Tile<D, S>::deposit_current() {
//
//  auto& yee = ::emf::Tile<D>::get_grids();
//
//  //std::cout << "Calling DAMPED J update\n";
//
//  //std::cout<<"dt:"<<yeeDt<<"  and vol:"<<yeeDx<<" .. " <<(yeeDx*yeeDy*yeeDz) <<"\n";
//  
//  float_m resistivity = 10.0;
//
//  yee.ex -= yee.jx / resistivity;
//  yee.ey -= yee.jy / resistivity;
//  yee.ez -= yee.jz / resistivity;
//
//}


/// Damp EM field into the reference field
//
// Combined directions:
//
// |S] == 1 : X
// |S] == 2 : Y
// |S] == 3 : Z
//
// -S is for left direction, +S for right
//
// NOTE: since S is known during compile time, compiler will optimize every other
// direction out but the correct one
template<
  std::size_t D, 
  int S
>
void emf::damping::Tile<D,S>::damp_fields()
{

  auto& yee = this->get_grids();

  float_m lambda2;

  bool left  = S < 0;
  bool right = S > 0;
    
  bool xdir  = abs(S) == 1;
  bool ydir  = abs(S) == 2;
  bool zdir  = abs(S) == 3;

  // fixed position
  const int halo = 1;

  int istr = 0;
  int ifin = emf::Tile<D>::mesh_lengths[0];

  int jstr = 0;
  int jfin = emf::Tile<D>::mesh_lengths[1];

  int kstr = 0;
  int kfin = emf::Tile<D>::mesh_lengths[2];

  // reset also halo regions
  if(D >= 1) istr -= halo;
  if(D >= 1) ifin += halo;

  if(D >= 2) jstr -= halo;
  if(D >= 2) jfin += halo;

  if(D >= 3) kstr -= halo;
  if(D >= 3) kfin += halo;


  float_m iglob, jglob, kglob, glob_coord;
  float_m t1, t2;

  bool inside_damping_zone = false;

  // loop over all values in mesh
  for(int k = kstr; k<kfin; k++) {
    kglob = (float_m)k + emf::Tile<D>::mins[2];

    for(int j = jstr; j<jfin; j++) {
      jglob = (float_m)j + emf::Tile<D>::mins[1];

      for(int i=istr; i<ifin; i++) {
        iglob = (float_m)i + emf::Tile<D>::mins[0];

        inside_damping_zone = false;

        // select active coordinate to reference against
        if(xdir) glob_coord = iglob;
        if(ydir) glob_coord = jglob;
        if(zdir) glob_coord = kglob;

        // left direction; 
        // fld1 < x < fld2 | <- incoming wave
        // <<< << <
        if(left  && ((glob_coord >= fld1) && (glob_coord < fld2))) {
          lambda2 = ksupp*pow(1.*(fld2-glob_coord)/(fld2-fld1), 3); 
          //lambda  = (1.0 - 0.5*lambda1)/(1.0 + 0.5*lambda1); //second order
          inside_damping_zone = true;

          t1 = exp(-lambda2);
          t2 = 1.0-exp(-lambda2);
        }

        // right direction; 
        // incoming wave -> | fld2 < x < fld1 
        // > >> >>>
        if(right && ((glob_coord <= fld1) && (glob_coord > fld2))) {
          lambda2 = ksupp*pow(1.*(fld1-glob_coord)/(fld1-fld2), 3);
          //lambda  = (1.0 - 0.5*lambda1)/(1.0 + 0.5*lambda1); //second order
          inside_damping_zone = true;

          t1 = 1.0-exp(-lambda2);
          t2 = exp(-lambda2);
        }

        if(inside_damping_zone){
            yee.ex(i,j,k) = t1*yee.ex(i,j,k) + t2*ex_ref(i,j,k);
            yee.ey(i,j,k) = t1*yee.ey(i,j,k) + t2*ey_ref(i,j,k);
            yee.ez(i,j,k) = t1*yee.ez(i,j,k) + t2*ez_ref(i,j,k);

            yee.bx(i,j,k) = t1*yee.bx(i,j,k) + t2*bx_ref(i,j,k);
            yee.by(i,j,k) = t1*yee.by(i,j,k) + t2*by_ref(i,j,k);
            yee.bz(i,j,k) = t1*yee.bz(i,j,k) + t2*bz_ref(i,j,k);
        }

        // fully damped to reference; behind the damping zone
        if ( (left &&  (glob_coord <= fld1)) 
          || (right && (glob_coord >  fld1)) ) {
          yee.ex(i,j,k) = ex_ref(i,j,k);
          yee.ey(i,j,k) = ey_ref(i,j,k);
          yee.ez(i,j,k) = ez_ref(i,j,k);

          yee.bx(i,j,k) = bx_ref(i,j,k);
          yee.by(i,j,k) = by_ref(i,j,k);
          yee.bz(i,j,k) = bz_ref(i,j,k);
        }
      }
    }
  }

}


//--------------------------------------------------
//--------------------------------------------------

// explicit template instantiation for supported directions
template class emf::damping::Tile<1,+1>;
template class emf::damping::Tile<1,-1>;

template class emf::damping::Tile<2,-1>;
template class emf::damping::Tile<2,+1>;
template class emf::damping::Tile<2,-2>;
template class emf::damping::Tile<2,+2>;
  
template class emf::damping::Tile<3,-1>;
template class emf::damping::Tile<3,+1>;
template class emf::damping::Tile<3,-2>;
template class emf::damping::Tile<3,+2>;
template class emf::damping::Tile<3,-3>;
template class emf::damping::Tile<3,+3>;
