/**
   @file    gridgas.cc
   @author  Luca Maccione, luca.maccione@desy.de
   @brief   Implementation of the TGridGas class. See the .h file
*/

#include "gridgas.h"

double TGridGas::Density(TVector3D aPosition) {

  // aPosition must be inside a shell around the grid
  double lLoc ;
  int lQx,lQy,lQz ;
  int lQx1,lQy1,lQz1 ;
  double a,b,c,d,e,f ;
  TVector3D lPosInStep = (aPosition-_fOrigin)/_fStepsize ;
  lQx = (int)floor(lPosInStep.x()) ;
  a=lPosInStep.x()-lQx ;
  d=1-a ;
  lQx1 = lQx+1 ;
  lQx = lQx % _fNx ;
  if (lQx<0) lQx += _fNx ;
  lQx1 = lQx1 % _fNx ;
  if (lQx1<0) lQx1 += _fNx ;
  lQy = (int)floor(lPosInStep.y()) ;
  b=lPosInStep.y()-lQy ;
  e=1-b ;
  lQy1 = lQy+1 ;
  lQy = lQy % _fNy ;
  if (lQy<0) lQy += _fNy ;
  lQy1 = lQy1 % _fNy ;
  if (lQy1<0) lQy1 += _fNy ;
  lQz = (int)floor(lPosInStep.z()) ;
  c=lPosInStep.z()-lQz ;
  f=1-c ;
  lQz1 = lQz+1 ;
  lQz = lQz % _fNz ;
  if (lQz<0) lQz += _fNz ;
  lQz1 = lQz1 % _fNz ;
  if (lQz1<0) lQz1 += _fNz ;
  
  int NyNz = _fNy*_fNz; 
  lLoc=d*(e*(f*_fpGas[lQx*NyNz+lQy*_fNz+lQz]+c*_fpGas[lQx*NyNz+lQy*_fNz+lQz1])
	   +b*(f*_fpGas[lQx*NyNz+lQy1*_fNz+lQz]+c*_fpGas[lQx*NyNz+lQy1*_fNz+lQz1]))+
    a*(e*(f*_fpGas[lQx1*NyNz+lQy*_fNz+lQz]+c*_fpGas[lQx1*NyNz+lQy*_fNz+lQz1])+
       b*(f*_fpGas[lQx1*NyNz+lQy1*_fNz+lQz]+c*_fpGas[lQx1*NyNz+lQy1*_fNz+lQz1])) ;
  
  if ( lLoc > this->Gasmax() ) 
    throw TCrpErr("Interpolated Gas field larger than Gasmax.") ;
  
  return lLoc;
}
