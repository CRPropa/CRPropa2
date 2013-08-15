/**
   @class TabulatedTALYSMeanFreePath_CMB

   @author  Nils Nierstenhoefer nierstenhoefer@physik.uni-wuppertal.de
 
   @brief Class to handle precalculated mean free path data for photo disintegration in the CMB (s. documentation of TabulatedTALYSY for details). 
*/

#ifndef _TabulatedTALYSMeanFreePath_CMB_H_
#define _TabulatedTALYSMeanFreePath_CMB_H_

#include <TabulatedTALYSMeanFreePath.h>

extern std::string TalysDirectory; // from sophiainteractions.cc

class TabulatedTALYSMeanFreePath_CMB : public TabulatedTALYSMeanFreePath{
public:
  static TabulatedTALYSMeanFreePath_CMB* GetInstance(){
    if(! fInstance_CMB){  
      fInstance_CMB = new  TabulatedTALYSMeanFreePath_CMB();
      fInstance_CMB->ReadData();
    }
    else { 
      //std::cout<<"TabulatedTALYSMeanFreePath_CMB already constructed."<<std::endl;
    }
    return fInstance_CMB; 
  }
protected:
  TabulatedTALYSMeanFreePath_CMB(){
    fPath<<(TalysDirectory+"CMB/");
    fPD_TabDataType=CMB_MFPTab;
  }
  ~TabulatedTALYSMeanFreePath_CMB();

private:
  static TabulatedTALYSMeanFreePath_CMB* fInstance_CMB;
  TabulatedTALYSMeanFreePath_CMB (const TabulatedTALYSMeanFreePath_CMB&);
  TabulatedTALYSMeanFreePath_CMB & operator = (const TabulatedTALYSMeanFreePath_CMB &);

};


#endif
