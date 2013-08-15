/**
   @class TabulatedTALYSMeanFreePath

   @author  Nils Nierstenhoefer nierstenhoefer@physik.uni-wuppertal.de
 
   @brief Class to handle precalculated mean free path data for photo disintegration (s. documentation of TabulatedTALYSY 
   for details). 

   No instance of this class is directly created - only used to derive TabulatedTALYSMeanFreePath_CMB and 
   TabulatedTALYSMeanFreePath_IRB via inheritance.  
*/

#ifndef _TabulatedTALYSMeanFreePath_H_
#define _TabulatedTALYSMeanFreePath_H_

#include "TabulatedTALYSY.h"
#include "TALYSMeanFreePath.h"
#include "PhotonBackground.h"
#include <TabulatedTALYSAveragedCrossSection.h>
#include <fstream>
#include<iostream>
#include <vector>
#include <algorithm>
#include <limits>

// This class is organized as a "factory" or "Singleton. That means you can only create one instance of this class which is accesible from all the programm. This is clearly useful for the cross-scetion data. It follows a short sketch which shows what the structure of a class with the needed which has this feature can look like:   
/*
class N
 {
 public:
    static N* instance()
    {
       if(! _instance)
       _instance = new N ();
       return _instance;
    }

    ~N();
    void xyz();

 private:
 static N* _instance;
    N();            // no instance of N can be created from the "outside" of N
                    // protected
    N( const N& );              //no new instance via copy constructor
    N & operator = (const N &); //no new instance via copy constructor
 };
*/

using namespace std;

class TabulatedTALYSMeanFreePath : public TabulatedTALYSY
{
 public:  
  /**< Defines which ExclY value is returned if x>FXmax.*/
  virtual double ExclYValueAboveFXMax(unsigned long int InputNucleusId,
				      unsigned long int ExclusiveChannelId){
    return(std::numeric_limits<double>::min());
  }

  /**< Defines which ExclY value is returned if x>FXmax.*/
  virtual double ExclSumYValueAboveFXMax(unsigned long int InputNucleusId){
    return(std::numeric_limits<double>::min());
  }
  
  virtual void CreateTables(std::vector<PhotonBackground*> PhotonBackgrounds_Vec,
			    std::vector<TabulatedTALYSY*> tTabulatedTALYSY_Vec,
			    PD_TabDataType tPD_TabDataType,
			    TVector3D Position,
			    double epsilon0, 
			    double epsilon_min, 
			    double epsilon_max,
			    std::string folder,
			    int NIntegratorCallsdouble,
			    double EMinMFPTableVal,
			    double EMaxMFPTableVal,
			    int NBinsMFPTableVal);

  virtual void ReduceTables(std::string folder,
			    double accuracy);
 protected:
  TabulatedTALYSMeanFreePath();
  void ReadData();
  std::string fDataPathExtension;
  std::ostringstream fPath;

 private:
  TabulatedTALYSMeanFreePath(const TabulatedTALYSMeanFreePath&);
  TabulatedTALYSMeanFreePath & operator = (const TabulatedTALYSMeanFreePath&);
  string fFolder;
  string fOption;
  void FillFiles(std::string folder,
		 unsigned long int  InitNucleusId,
		 std::vector< std::vector<long unsigned int> > * ExclusiveChannelImprtortant,	  
		 unsigned long int* LineCounterOutputNuclei,
		 unsigned long int* LineCounterOutExclChannel,
		 unsigned long int* LineCounterSumExclChannel,
		 unsigned long int* ExclCrossSecLineCounter,
		 ofstream* ExclusiveCrossSecData,
		 ofstream* ExclusiveCrossSecDataId,
		 ofstream* TALYSInitNucleusId,
		 ofstream* ExclusiveSumCrossSecData);
};
#endif

