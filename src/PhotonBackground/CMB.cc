
#include "CMB.h"

#include <cmath>

//For CMB x,y and z aren't needed!
//h_bar=c=1


CMB::CMB() {
  //fEpsiolonZero=k*T  
  fEpsilonZero= 8.617e-14*2.725;         
}

CMB::~CMB(){}

double CMB::GetPhotonDensity(double x, 
			     double y, 
			     double z, 
			     double redshift,
			     double PhotonEnergy){


  return(pow((1.+redshift), 2)* 
	 PhotonEnergy * PhotonEnergy / ( exp( PhotonEnergy / (fEpsilonZero*(1.+redshift))  ) 
					 - 1. ) 
	 /  
	 M_PI / M_PI 
	  );
}


