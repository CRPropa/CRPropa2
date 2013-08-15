/**
   @class TALYSMeanFreePathTabulated

   @author  Nils Nierstenhoefer nierstenhoefer@physik.uni-wuppertal.de
 
   @brief Class to directly read the mean free path for photo disintegration from the mean free path data tables 
   TabulatedTALYSMeanFreePath_CMB and TabulatedTALYSMeanFreePath_IRB (s. documentation of TALYSMeanFreePath for details). 
*/

#ifndef _TALYSMeanFreePathTabulated_H_
#define _TALYSMeanFreePathTabulated_H_

#include <cmath>
#include <vector>
#include <TabulatedTALYSCrossSections.h>
#include <iostream>
#include <iomanip>
#include <CMB.h>
#include <TALYSMeanFreePath.h> 
#include <TabulatedTALYSMeanFreePath.h>
#include <vector>



class TALYSMeanFreePathTabulated : public TALYSMeanFreePath{
/* 
typedef double (TALYSMeanFreePathTabulated::*ClassMemPtr) (double,
						   unsigned long int,
						   unsigned long int,
						   double) ;
*/

public:
  
  
  typedef double (TALYSMeanFreePathTabulated::*ClassMemPtr) (double,
							     unsigned long int,
							     unsigned long int,
							     double,
							     double,
							     double,
							     int);
  
  //  typedef double (TabulatedTALYSCrossSection::*TabulatedTALYSCrossSectionMemPntr) (
  //						    unsigned long int,
  //						    unsigned long int,
  //						    int);

  TALYSMeanFreePathTabulated();
  ~TALYSMeanFreePathTabulated();
  /*virtual double Total_rate(unsigned long int InputNucleusId,
			    unsigned long int OutPutNucleusId,
			    double NucleusMass,          //GeV/c^2
			    double EnergyOfNucleusGeV,   //GeV
			    double z,
			    int TotalOrExclusive
			    );*/
  
  virtual double Total_rate(std::vector<PhotonBackground*> PhotonBkgrd,
			    std::vector<TabulatedTALYSY*> tTabulatedTALYSY_Vec,
			    PD_TabDataType tPD_TabDataType,
			    TVector3D Position,
			    unsigned long int InputNucleusId,
			    unsigned long int OutPutNucleusId,
			    double NucleusMass,
			    double EnergyOfNucleusGeV,
			    double z,
			    int TotalOrExclusive);
  
  virtual double Total_rate_log(std::vector<PhotonBackground*> PhotonBkgrd,
				std::vector<TabulatedTALYSY*> tTabulatedTALYSY_Vec,
				PD_TabDataType tPD_TabDataType,
				TVector3D Position,
				unsigned long int InputNucleusId,
				unsigned long int OutPutNucleusId,
				double NucleusMass,
				double EnergyOfNucleusGeV,
				double z,
				int TotalOrExclusive){

    //In case of the mean free path tables Total_rate() and Total_rate_log() should return the same results
    return(Total_rate(PhotonBkgrd,
		      tTabulatedTALYSY_Vec,
		      tPD_TabDataType,
		      Position,
		      InputNucleusId,
		      OutPutNucleusId,
		      NucleusMass,
		      EnergyOfNucleusGeV,
		      z,
		      TotalOrExclusive));
  }

  virtual double functs(double s, 
			TVector3D Position,
			unsigned long int InputNucleusId,
			unsigned long int OutputId,
			double NucleusMass,
			double EnergyOfNucleusGeV,
			double z,
			int TotalOrExclusive);

  /*  virtual double functionToTestIntegrator(double s, 
		unsigned long int InputNucleusId,
		unsigned long int OutputId,
		double NucleusMass,
		double EnergyOfNucleusGeV,
		double z,
		int TotalOrExclusive);
  */
 protected:

 private:

};

#endif
