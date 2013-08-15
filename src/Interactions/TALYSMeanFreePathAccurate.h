/**
   @class TALYSMeanFreePathAccurate

   @author  Nils Nierstenhoefer nierstenhoefer@physik.uni-wuppertal.de
 
   @brief Class to directly calculate the mean free path for photo disintegration from the cross section data TabulatedTALYSCrossSection 
   (s. documentation of TALYSMeanFreePath for details). 

   Mainly used to calculate the "averged cross section" data (c.f. TabulatedTALYSAveragedCrossSection).

   Note, the term "accurate" in this classname is missleading and should be renamed to e.g. TALYSMeanFreePathFromCrossSection. 
*/

#ifndef _TALYSMeanFreePathAccurate_H_
#define _TALYSMeanFreePathAccurate_H_

#include <cmath>
#include <vector>
#include <TabulatedTALYSCrossSections.h>
#include <iostream>
#include <iomanip>
#include <CMB.h>
#include <TALYSMeanFreePath.h> 

class TALYSMeanFreePathAccurate : public TALYSMeanFreePath{
/* 
typedef double (TALYSMeanFreePathAccurate::*ClassMemPtr) (double,
						   unsigned long int,
						   unsigned long int,
						   double) ;
*/

public:
  
  
  typedef double (TALYSMeanFreePathAccurate::*ClassMemPtr) (double,
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

  TALYSMeanFreePathAccurate();
  ~TALYSMeanFreePathAccurate();
  
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

  virtual double functs_int(double UpperBorder,
		    unsigned long int InputNucleusId,
		    unsigned long int OutPutNucleusId,
		    double NucleusMass,
		    double EnergyOfNucleusGeV,
		    double z,
		    int TotalOrExclusive);
  virtual double functs_int_times_PhtnBckGrnd(double t,
				       unsigned long int InputNucleusId,
				       unsigned long int OutPutNucleusId,
				       double NucleusMass,
				       double EnergyOfNucleusGeV,
				       double z,
				       int TotalOrExclusive);
  virtual double Total_rate_cmb(unsigned long int InputNucleusId,
				unsigned long int OutPutNucleusId,
				double NucleusMass,
				double EnergyOfNucleusGeV,
				double z,
				int TotalOrExclusive);
  //double gauss (ClassMemPtr function , double A, double B);
  //double PHOTD(double EPS,double TBB);
  //double PHOTD(double t);
  
 protected:
  virtual double gauss (ClassMemPtr function,
			double A,
			double B,
			unsigned long int InputNucleusId,
			unsigned long int OutPutNucleusId,
			double NucleusMass,
			double EnergyOfNucleusGeV,
			double z,
			int TotalOrExclusive
			);
  //Variables
  //Create an instance of TabulatedTALYSCrossSection
  //TabulatedTALYSCrossSection fTabulatedTALYSCrossSection;
  */
};

#endif
