/**
   @class TabulatedTALYSCrossSection

   @author  Nils Nierstenhoefer nierstenhoefer@physik.uni-wuppertal.de

   @brief Class to handle photo disintegration cross section data (s. documentation of TabulatedTALYSY for details).  
*/

#ifndef _TabulatedTALYSCrossSection_H_
#define _TabulatedTALYSCrossSection_H_

#include "TabulatedTALYSY.h"

#include <fstream>
#include<iostream>
#include <vector>
#include <algorithm>
#include <typedclass.h>

#include <TGraph.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TLegend.h>


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

class TabulatedTALYSCrossSection : public TabulatedTALYSY 
{
 public:


  static TabulatedTALYSCrossSection* GetInstance(){ 
    if(! finstance) 
      finstance = new TabulatedTALYSCrossSection(); 
    return finstance; 
  } 
  ~TabulatedTALYSCrossSection();

  std::string PolymorphieTester() {return("TabulatedTALYSCrossSection");}
 
  virtual double GetTabExclSummedYValue(unsigned long int InputNucleusId,
					int EnergyBin){
    throw TCrpErr("Error: TabulatedTALYSCrossSection::GetTabExclSummedYValue(): not implemented in this context.");
  }
  
  virtual double GetExclusiveSummedY(unsigned long int InputNucleusId,
				     double Energy){
    throw TCrpErr("Error: TabulatedTALYSCrossSection::GetTabExclSummedYValue(): not implemented in this context.");
  }

 protected:

 private:
  static TabulatedTALYSCrossSection* finstance;
  TabulatedTALYSCrossSection();
  TabulatedTALYSCrossSection(const TabulatedTALYSCrossSection&);
  TabulatedTALYSCrossSection & operator = (const TabulatedTALYSCrossSection&);
};
#endif

