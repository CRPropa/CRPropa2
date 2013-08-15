
#include "Opt.h"
#include <cmath>

Opt::Opt() {
  //fEpsilonZero= 0;
#ifdef DEBUG_OUTPUT
  //std::cout<<"WARNING: Optical background is not implemented!"<<std::endl;
#endif
}

Opt::~Opt(){}

double Opt::GetPhotonDensity(double x, 
			     double y, 
			     double z, 
			     double redshift,
			     double PhotonEnergy){

  return(0.);
}


