
/**
   @file    particle.cc
   @author  Eric Armengaud, armengaud@in2p3.fr
   @brief   Generic class for particles
*/

#include "particle.h"
#include <algorithm>
#include <vector>

QUEUE<TParticle*> TParticle::fParticleQueue ;


void TParticle::ComputeRedshift() {

  if (_fpUniverse->Type() != UNIVERSE_ENV1D ){
    _fRedshift = 0 ;
  } else {
    unsigned int i0 = 0 ;
      //std::vector<double> da = _fpUniverse->DistanceArray();
      //i0 = (upper_bound(da.begin(),da.end(),_fPosition.x())-da.begin()-1);
      //cout << "Begin " << endl;
      //cout << "i0 new = " << i0 << endl;
      i0=0;
      while (_fpUniverse->DistanceArray(i0+1) <= _fPosition.x() ) i0 += 1 ;
      //cout << "i0 old = " << i0 << endl;
    _fRedshift = _fpUniverse->RedshiftArray(i0) + (_fpUniverse->RedshiftArray(i0+1)-_fpUniverse->RedshiftArray(i0))
      *(_fPosition.x()-_fpUniverse->DistanceArray(i0))
      /(_fpUniverse->DistanceArray(i0+1)-_fpUniverse->DistanceArray(i0)) ;
  if (_fPosition.x() <= 0) _fRedshift = 0 ;
  }

}
/*
void TParticle::ComputeRedshift(const TUniverse* aUniv) {

  if (aUniv->Type() != UNIVERSE_ENV1D) {
    _fRedshift = 0 ;
  } else {
    unsigned int i0 = 0 ;
    while (aUniv->DistanceArray(i0+1) <= _fPosition.x() ) i0 += 1 ;
    _fRedshift = aUniv->RedshiftArray(i0) + (aUniv->RedshiftArray(i0+1)-aUniv->RedshiftArray(i0))
      *(_fPosition.x()-aUniv->DistanceArray(i0))
      /(aUniv->DistanceArray(i0+1)-aUniv->DistanceArray(i0)) ;
  if (_fPosition.x() <= 0) _fRedshift = 0 ;
  }

}
*/
