/**
   @class TALYSMeanFreePathAvrgd

   @author  Nils Nierstenhoefer nierstenhoefer@physik.uni-wuppertal.de
 
   @brief Class to directly calculate the mean free path for photo disintegration from the "averaged cross section" data 
   TabulatedTALYSAveragedCrossSection (s. documentation of TALYSMeanFreePath for details). 
*/

#ifndef _TALYSMeanFreePathAvrgd_H_
#define _TALYSMeanFreePathAvrgd_H_

#include <cmath>
#include <vector>
#include <TabulatedTALYSCrossSections.h>
#include <iostream>
#include <iomanip>
#include <CMB.h>
#include <TALYSMeanFreePath.h>

class TALYSMeanFreePathAvrgd : public TALYSMeanFreePath{
/* 
typedef double (TALYSMeanFreePathAvrgd::*ClassMemPtr) (double,
						   unsigned long int,
						   unsigned long int,
						   double) ;
*/

public:
  
  
  typedef double (TALYSMeanFreePathAvrgd::*ClassMemPtr) (double,
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
  TALYSMeanFreePathAvrgd(double epsilon_zero, double epsilon_min, double epsilon_max);
  ~TALYSMeanFreePathAvrgd();
  
  /*
  virtual double Total_rate(unsigned long int InputNucleusId,
  unsigned long int OutPutNucleusId,
			    double NucleusMass,
			    double EnergyOfNucleusGeV,
			    double z,
			    int TotalOrExclusive){
			    std::cerr<<"Error: Choose another overloaded implementation of Total_rate(...) which calculates MFP on the fly."<<std::endl;
    } 
  */
  
   std::string PolymorphieTester() {return ("TALYSMeanFreePathAvrgd");}
  
   virtual double functs_int(double UpperBorder,
			     TVector3D Position,
			     unsigned long int InputNucleusId,
			     unsigned long int OutPutNucleusId,
			     double NucleusMass,
			     double EnergyOfNucleusGeV,
			     double z,
			     int TotalOrExclusive);

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

 private:
  TALYSMeanFreePathAvrgd();
  
};

#endif
