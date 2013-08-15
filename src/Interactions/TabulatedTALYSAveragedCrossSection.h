/**
   @class TabulatedTALYSAveragedCrossSection

   @author  Nils Nierstenhoefer nierstenhoefer@physik.uni-wuppertal.de

   @brief Class to handle "averaged cross section" data for photo disintegration (s. documentation of TabulatedTALYSY for details and the corresponding  definition).  
*/

#ifndef _TabulatedTALYSAveragedCrossSection_H_
#define _TabulatedTALYSAveragedCrossSection_H_

#include "TabulatedTALYSY.h"
#include "TabulatedTALYSCrossSections.h"

#include <fstream>
#include<iostream>
#include <vector>
#include <algorithm>

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
                    // protected, wenn man von der Klasse noch erben möchte
    N( const N& );              //no new instance via copy constructor
    N & operator = (const N &); //no new instance via copy constructor
 };
*/

using namespace std;

class TabulatedTALYSAveragedCrossSection : public TabulatedTALYSY 
{
 public:

  static TabulatedTALYSAveragedCrossSection* GetInstance(){
     if(! finstance)
       finstance = new TabulatedTALYSAveragedCrossSection();
       return finstance;
  } 
  ~TabulatedTALYSAveragedCrossSection();
  
  virtual std::string PolymorphieTester() {return("TabulatedTALYSAveragedCrossSection");}

  void CreateTables(std::string folder);
  void ReduceTables(std::string folder,
		    double accuracy);

  /**< Defines which ExclY value is returned if x>FXmax.*/
  double ExclYValueAboveFXMax(unsigned long int InputNucleusId,
					  unsigned long int ExclusiveChannelId){
    return(GetTabExclChannelYValue(InputNucleusId,
				   ExclusiveChannelId,
				   GetNXBins()-1));
  }
  
  /**< Defines which ExclY value is returned if x>FXmax.*/
  double ExclSumYValueAboveFXMax(unsigned long int InputNucleusId){
    return(GetTabExclSummedYValue(InputNucleusId,
				  GetNXBins()-1));
  }
  
 protected:

 private:
  TabulatedTALYSAveragedCrossSection();
  static TabulatedTALYSAveragedCrossSection* finstance;
  TabulatedTALYSAveragedCrossSection(const TabulatedTALYSAveragedCrossSection&);
  TabulatedTALYSAveragedCrossSection & operator = (const TabulatedTALYSAveragedCrossSection&);
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

