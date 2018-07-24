#include "damping_fields.h"



using std::min;
using std::max;
using std::exp;


template<int S>
fields::PlasmaCellDamped<S>::PlasmaCellDamped(
    size_t i, size_t j,
    int o,
    size_t NxG, size_t NyG,
    size_t NxMesh, size_t NyMesh, size_t NzMesh
    ) : 
  corgi::Cell(i, j, o, NxG, NyG),
  PlasmaCell(i, j, o, NxG, NyG, NxMesh, NyMesh, NzMesh),
  ex_ref(NxMesh, NyMesh, NzMesh),
  ey_ref(NxMesh, NyMesh, NzMesh),
  ez_ref(NxMesh, NyMesh, NzMesh),

  bx_ref(NxMesh, NyMesh, NzMesh),
  by_ref(NxMesh, NyMesh, NzMesh),
  bz_ref(NxMesh, NyMesh, NzMesh)
{ }



/*
void fields::PlasmaCellDamped::pushE() {

  // this->pushE1d();
  this->pushE2d_damped();
  // this->pushE3d();
}
*/


/*
/// 2D E pusher
void fields::PlasmaCellDamped::pushE2d_damped() {

  fields::YeeLattice& mesh = getYee();

  //std::cout << "Calling DAMPED E update\n";
  Realf C = 1.0 * cfl;
 
  int k = 0;
  for(int j=0; j<(int)NyMesh; j++) {
    for(int i=0; i<(int)NxMesh; i++) {

      // Ex
      mesh.ex(i,j,k) += 
        + C*(-mesh.bz(i,j-1,k  ) + mesh.bz(i,j,k));

      // Ey
      mesh.ey(i,j,k) += 
        + C*( mesh.bz(i-1,j, k  ) - mesh.bz(i,j,k));

      // Ez
      mesh.ez(i,j,k) += 
        + C*( mesh.bx(i,  j-1, k) - mesh.bx(i,j,k))
        + C*(-mesh.by(i-1,j,   k) + mesh.by(i,j,k));

    }
  }
}
*/


// deposit current with some resistivity
template<int S>
void fields::PlasmaCellDamped<S>::depositCurrent() {
  fields::YeeLattice& mesh = getYee();

  //std::cout << "Calling DAMPED J update\n";

  //std::cout<<"dt:"<<yeeDt<<"  and vol:"<<yeeDx<<" .. " <<(yeeDx*yeeDy*yeeDz) <<"\n";
  
  Realf resistivity = 10.0;

  mesh.ex -= mesh.jx / resistivity;
  mesh.ey -= mesh.jy / resistivity;
  mesh.ez -= mesh.jz / resistivity;

}


/// Damp EM fields into the reference field
//
// X-dir
template<int S>
template<int R>
typename std::enable_if<R==-1 | R==+1>::type 
fields::PlasmaCellDamped<S>::dampFields()
{
  auto& yee = getYee();

  //Realf lambda, lambda1;
  Realf lambda2;
  Realf ksupp = 10; // damping parameter


  // moving injector position
  //int ntimes = 10; // filter width
  //int iback = 10 + max(50, ntimes);
  //int istr = max(0,  (int)(xinject-iback) );
  //int ifin = min(mx, (int)(xinject+iback) );
    
  // fixed position
  int istr = 0;
  int ifin = NxMesh;

  int jstr = 0;
  int jfin = NyMesh;


  int k = 0; // XXX: 2D hack
  Realf iglob, jglob, kglob;

  kglob = (Realf)k + mins[2];
  for(int i = istr; i<ifin; i++) {
    iglob = (Realf)i + mins[0];

    // damping region
    if ( ((S<0) && ((iglob >= fld1) && (iglob < fld2))) || 
         ((S>0) && ((iglob <= fld1) && (iglob > fld2))) )  
    {
      //lambda1 = kappa*pow( (fld2-jglob)/(fld2 - fld1), 3);
      //lambda  = (1.0 - 0.5*lambda1)/(1.0 + 0.5*lambda1); //second order

      if (S < 0) lambda2 = ksupp*pow(1.*(fld2-iglob)/(fld2-fld1), 3); // left dir
      if (S > 0) lambda2 = ksupp*pow(1.*(fld1-iglob)/(fld1-fld2), 3); // right dir

      for(int j=jstr; j<jfin; j++) {
        jglob = (Realf)j + mins[1];

        yee.ex(i,j,k) = exp(-lambda2)*yee.ex(i,j,k) + (1.0 - exp(-lambda2))*ex_ref(i,j,k);
        yee.ey(i,j,k) = exp(-lambda2)*yee.ey(i,j,k) + (1.0 - exp(-lambda2))*ey_ref(i,j,k);
        yee.ez(i,j,k) = exp(-lambda2)*yee.ez(i,j,k) + (1.0 - exp(-lambda2))*ez_ref(i,j,k);

        yee.bx(i,j,k) = exp(-lambda2)*yee.bx(i,j,k) + (1.0 - exp(-lambda2))*bx_ref(i,j,k);
        yee.by(i,j,k) = exp(-lambda2)*yee.by(i,j,k) + (1.0 - exp(-lambda2))*by_ref(i,j,k);
        yee.bz(i,j,k) = exp(-lambda2)*yee.bz(i,j,k) + (1.0 - exp(-lambda2))*bz_ref(i,j,k);

      }
    }

    // already damped region
    // left || right
    if( ((S < 0) && (iglob < fld1)) || ((S > 0) && (iglob > fld1)) ) {
      for(int j=jstr; j<jfin; j++) {
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


// Y-dir
template<int S>
template<int R>
typename std::enable_if<R==-2 | R==+2>::type 
fields::PlasmaCellDamped<S>::dampFields()
{
  auto& yee = getYee();

  //Realf lambda, lambda1;
  Realf lambda2;
  Realf ksupp = 10; // damping parameter


  // moving injector position
  //int ntimes = 10; // filter width
  //int iback = 10 + max(50, ntimes);
  //int istr = max(0,  (int)(xinject-iback) );
  //int ifin = min(mx, (int)(xinject+iback) );
    
  // fixed position
  int istr = 0;
  int ifin = NxMesh;


  int k = 0; // XXX: 2D hack
  Realf jglob;
  for(int j = 0; j<(int)NyMesh; j++) {
    jglob = (Realf)j + mins[1];

    // damping region
    if ( ((S<0) && ((jglob >= fld1) && (jglob < fld2))) || 
         ((S>0) && ((jglob <= fld1) && (jglob > fld2))) )  
    {
      //lambda1 = kappa*pow( (fld2-jglob)/(fld2 - fld1), 3);
      //lambda  = (1.0 - 0.5*lambda1)/(1.0 + 0.5*lambda1); //second order

      if (S < 0) lambda2 = ksupp*pow(1.*(fld2-jglob)/(fld2-fld1), 3);
      if (S > 0) lambda2 = ksupp*pow(1.*(fld1-jglob)/(fld1-fld2), 3);

      for(int i=istr; i<ifin; i++) {
        yee.ex(i,j,k) = exp(-lambda2)*yee.ex(i,j,k) + (1.0 - exp(-lambda2))*ex_ref(i,j,k);
        yee.ey(i,j,k) = exp(-lambda2)*yee.ey(i,j,k) + (1.0 - exp(-lambda2))*ey_ref(i,j,k);
        yee.ez(i,j,k) = exp(-lambda2)*yee.ez(i,j,k) + (1.0 - exp(-lambda2))*ez_ref(i,j,k);

        yee.bx(i,j,k) = exp(-lambda2)*yee.bx(i,j,k) + (1.0 - exp(-lambda2))*bx_ref(i,j,k);
        yee.by(i,j,k) = exp(-lambda2)*yee.by(i,j,k) + (1.0 - exp(-lambda2))*by_ref(i,j,k);
        yee.bz(i,j,k) = exp(-lambda2)*yee.bz(i,j,k) + (1.0 - exp(-lambda2))*bz_ref(i,j,k);

      }
    }

    // already damped region
    if( ((S < 0) && (jglob < fld1)) || ((S > 0) && (jglob > fld1)) ) {
      for(int i=istr; i<ifin; i++) {
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


//--------------------------------------------------
//--------------------------------------------------

// explicit template instantiation for supported directions
  
template class fields::PlasmaCellDamped<-1>;
template class fields::PlasmaCellDamped<+1>;
template class fields::PlasmaCellDamped<-2>;
template class fields::PlasmaCellDamped<+2>;
  
//template void  fields::PlasmaCellDamped<-1>::dampFields();
//template void  fields::PlasmaCellDamped<+1>::dampFields();
//template void  fields::PlasmaCellDamped<-2>::dampFields();
//template void  fields::PlasmaCellDamped<+2>::dampFields();

