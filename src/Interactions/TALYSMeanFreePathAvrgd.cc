#include "TALYSMeanFreePathAvrgd.h"
#include "TabulatedTALYSAveragedCrossSection.h"

TALYSMeanFreePathAvrgd::TALYSMeanFreePathAvrgd(double epsilon_zero, double epsilon_min, double epsilon_max)
: TALYSMeanFreePath(epsilon_zero, epsilon_min, epsilon_max){
}

TALYSMeanFreePathAvrgd::~TALYSMeanFreePathAvrgd(){}

//integrand for Total_rate_cmb / inner of the double integral from table
double TALYSMeanFreePathAvrgd::functs_int(double UpperBorder,
					  TVector3D Position,
					  unsigned long int InputNucleusId,
					  unsigned long int OutPutNucleusId,
					  double NucleusMass,
					  double EnergyOfNucleusGeV,
					  double z,
					  int TotalOrExclusive){
  if (UpperBorder < 0.) return(0.); 
  
  if(TotalOrExclusive==1){
    return(
	   TabulatedTALYSAveragedCrossSection::GetInstance()->
	   GetExclusiveChannelY(InputNucleusId,
				OutPutNucleusId,
				UpperBorder*1000));
    
  }else if(TotalOrExclusive==2){
    return(
	   TabulatedTALYSAveragedCrossSection::GetInstance()->
	   GetExclusiveSummedY(InputNucleusId,
			       UpperBorder*1000));
    
  }else{
    throw TCrpErr("Error in TALYSMeanFreePathAvrgd::functs(): TotalOrExclusive out of bound."); 
  }
}


//returns (s-pm^2)*sigma_Ngamma which corresponds to epsprime*sigma(epsprime) 
double TALYSMeanFreePathAvrgd::functs(double s, 
				      TVector3D Position,
				      unsigned long int InputNucleusId,
				      unsigned long int OutputId,
				      double NucleusMass,
				      double EnergyOfNucleusGeV,
				      double z,
				      int TotalOrExclusive){
  
  std::cout<<"TALYSMeanFreePathAvrgd::functs() called. This function is not needed here because we use the TabultedTALYSAveragedCrossSections."<<std::endl; 
  exit(-1);
}

