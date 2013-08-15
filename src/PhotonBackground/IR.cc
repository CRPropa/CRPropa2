#include "IR.h"
#include <iostream>
#include <cmath>
#include <vector>


//Give PhotonEnergy in Gev. Returns Density in GeV^2
double IR::GetPhotonDensity(double x, 
			    double y, 
			    double z, 
			    double redshift,
			    double PhotonEnergy){


  //PRIMACK
  std::vector<double> XData;
  XData.push_back(-1.); 
  XData.push_back(-0.75); 
  XData.push_back(-0.5); 
  XData.push_back(-0.25); 
  XData.push_back(0.);
  XData.push_back(0.25); 
  XData.push_back(0.5); 
  XData.push_back(0.75); 
  XData.push_back(1.); 
  XData.push_back(1.25); 
  XData.push_back(1.5); 
  XData.push_back(1.75); 
  XData.push_back(2.); 
  XData.push_back(2.25); 
  XData.push_back(2.5);

  std::vector<double> YData;
  YData.push_back(0.8); 
  YData.push_back(1.1); 
  YData.push_back(1.15); 
  YData.push_back(1.2); 
  YData.push_back(1.3);
  YData.push_back(1.2); 
  YData.push_back(1.05); 
  YData.push_back(0.7); 
  YData.push_back(0.4); 
  YData.push_back(0.3); 
  YData.push_back(0.5); 
  YData.push_back(0.8); 
  YData.push_back(1.1); 
  YData.push_back(1.3); 
  YData.push_back(1.);

  double REDSHIFT, Z, ZMAX_IR=_fZmax;
  
  //Redshift evolution as for CMB. (added Dec.'05)
  //conversion from nW/cm^3/sr to eV/cm^3
  //        FLUX_CONVERSION = 2.9979e10/4./dacos(-1.)*1.602e-19*1.e9*1.e4
  double FLUX_CONVERSION = 3.82182e3;
  
  //conversion GeV -> eV
  PhotonEnergy=PhotonEnergy*1.e9;

  //c conversion from eV to micrometers
  double X = 1.2398*(1.+redshift)/PhotonEnergy;
  double RESULT=0.;
  if (X>500.){
    RESULT=0.;
  }   
  else if(log10(X)>=XData[15-1]){ 
    RESULT = (YData[15-1] - YData[14-1])/(XData[15-1] - XData[14-1])*(log10(X) - XData[14-1]) + YData[14-1];
    RESULT = pow(10.,RESULT);
  }
  
  if (log10(X) <= XData[1-1]) RESULT=0.;
  
  if ((log10(X) < XData[15-1]) && (log10(X) > XData[1-1])){
    int INDEX = 2;
    do 
      INDEX = INDEX+1;
    while (XData[INDEX-1] < log10(X));
    
    RESULT = (YData[INDEX-1] - YData[INDEX-1-1])/ (XData[INDEX-1] - XData[INDEX-1-1])*(log10(X) - XData[INDEX-1-1])+  YData[INDEX-1-1];
    
    RESULT = pow(10., RESULT);  
  }
  
  if (Z>ZMAX_IR) {
    return(0.);
  }
  
  return(RESULT*pow((1.+redshift), 2)/pow((PhotonEnergy/(1.+redshift)), 2)/FLUX_CONVERSION/pow(5.068,3)/1.e30);
}
  
  
