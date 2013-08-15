#include "TALYSMeanFreePathAccurate.h"

TALYSMeanFreePathAccurate::TALYSMeanFreePathAccurate(){}

TALYSMeanFreePathAccurate::~TALYSMeanFreePathAccurate(){}

//returns (s-pm^2)*sigma_Ngamma which corresponds to epsprime*sigma(epsprime) 
double TALYSMeanFreePathAccurate::functs(double s, 
					 TVector3D Position,
					 unsigned long int InputNucleusId,
					 unsigned long int OutputId,
					 double NucleusMass,
					 double EnergyOfNucleusGeV,
					 double z,
					 int TotalOrExclusive){  
  double crossection;  
  double factor = s;
  double epsprime = factor;
  double sigma_pg;
  
  if(TotalOrExclusive==1){
    sigma_pg = 
      TabulatedTALYSCrossSection::GetInstance()->GetExclusiveChannelY(InputNucleusId,
								      OutputId,
								      epsprime*1000.);
  }else if(TotalOrExclusive==2){
    throw TCrpErr("Error in TALYSMeanFreePathAccurate::functs(): TotalOrExclusive==2 is not implemented in this context."); 
  }else{
    throw TCrpErr("Error in TALYSMeanFreePathAccurate::functs(): TotalOrExclusive out of bound."); 
  }
  
  return(epsprime*sigma_pg);
}
