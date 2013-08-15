/**
   @class TabulatedTALYSMeanFreePath_IRB

   @author  Nils Nierstenhoefer nierstenhoefer@physik.uni-wuppertal.de
 
   @brief Class to handle precalculated mean free path data for photo disintegration in the IRB (s. documentation of TabulatedTALYSY for details).
*/

#ifndef _TabulatedTALYSMeanFreePath_IRB_H_
#define _TabulatedTALYSMeanFreePath_IRB_H_

#include <TabulatedTALYSMeanFreePath.h>

extern std::string TalysDirectory; // from sophiainteractions.cc

class TabulatedTALYSMeanFreePath_IRB : public TabulatedTALYSMeanFreePath{
public:
  static TabulatedTALYSMeanFreePath_IRB* GetInstance(){
    if(! fInstance_IRB){  
      fInstance_IRB = new  TabulatedTALYSMeanFreePath_IRB();
      fInstance_IRB->ReadData();
    }
    else { 
      //std::cout<<"TabulatedTALYSMeanFreePath_IRB already constructed."<<std::endl;
    }
    return fInstance_IRB; 
  }
protected:
  TabulatedTALYSMeanFreePath_IRB(){
    fPath<<(TalysDirectory+"IRB/");
    fPD_TabDataType=IRB_MFPTab;
  }
  ~TabulatedTALYSMeanFreePath_IRB();

private:
  static TabulatedTALYSMeanFreePath_IRB* fInstance_IRB;
  TabulatedTALYSMeanFreePath_IRB (const TabulatedTALYSMeanFreePath_IRB&);
  TabulatedTALYSMeanFreePath_IRB & operator = (const TabulatedTALYSMeanFreePath_IRB &);
};

#endif
