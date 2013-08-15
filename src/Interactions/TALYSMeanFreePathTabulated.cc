#include "TALYSMeanFreePathTabulated.h"
#include "TabulatedTALYSMeanFreePath_CMB.h"
#include "TabulatedTALYSMeanFreePath_IRB.h"


TALYSMeanFreePathTabulated::TALYSMeanFreePathTabulated(){
}

TALYSMeanFreePathTabulated::~TALYSMeanFreePathTabulated(){}

//RETURNS total interaction rate with CMB in Mpc^-1
double TALYSMeanFreePathTabulated::Total_rate(std::vector<PhotonBackground*> PhotonBkgrd,
					      std::vector<TabulatedTALYSY*> tTabulatedTALYSY_Vec,
					      PD_TabDataType tPD_TabDataType,
					      TVector3D Position,
					      unsigned long int InputNucleusId,
					      unsigned long int OutPutNucleusId,
					      double NucleusMass,          //GeV/c^2
					      double EnergyOfNucleusGeV,   //GeV
					      double z,
					      int TotalOrExclusive
					      ){
  
  //sanity check . . . there should be no photonfields given in that case as we use the mfp tables.
  if(PhotonBkgrd.size()!=0){
    throw TCrpErr("TALYSMeanFreePathTabulated::Total_rate(...): Photonfields were given but can not be handeld.");
  }
  if(tTabulatedTALYSY_Vec.size()==0 || tTabulatedTALYSY_Vec.size()>2) throw TCrpErr("Total_rate(...): tTabulatedTALYSY_Vec.size()=0 or >2!"); 
  if(tPD_TabDataType == AveragedXSTab 
     || 
     tPD_TabDataType == XSTab) throw TCrpErr("Total_rate(...): cross sections or averaged cross sections data given but MFP data nedded."); 
  
  //find component in vector tTabulatedTALYSY_Vec which holds the demanded data of type tPD_TabDataType
  int tTabulatedTALYSY_Vec_Component = 0;
  int NComponents = tTabulatedTALYSY_Vec.size(); 
  if(NComponents>0){
    for(int a=0; a<NComponents; a++){
      tTabulatedTALYSY_Vec_Component=a;
      if (tTabulatedTALYSY_Vec[a]->GetPD_TabDataType()==tPD_TabDataType) break;
    }
  }

  double  xmp = NucleusMass;   

  double LogGamma = log10(EnergyOfNucleusGeV/NucleusMass);

  if(z>0.){
    std::cout<<"Warning in TALYSMeanFreePath::Total_rate_cmb(...): z=0 as you chose to use the tabulated MFP."<<std::endl;
    std::cout<<"z="<<z<<std::endl;
  }       
    
    
    double MFP = 0.0;
    switch (TotalOrExclusive) {
        case 1:
            MFP = tTabulatedTALYSY_Vec[tTabulatedTALYSY_Vec_Component]->GetExclusiveChannelY(InputNucleusId,
                                                                                             OutPutNucleusId,
                                                                                             LogGamma);
            break;
        case 2:
            MFP = tTabulatedTALYSY_Vec[tTabulatedTALYSY_Vec_Component]->GetExclusiveSummedY(InputNucleusId,
                                                                                            LogGamma);
            break;
        default:
            throw TCrpErr("Error in TALYSMeanFreePathTabulated::functs(): TotalOrExclusive out of bound.");
            break;
    }
    //check if inf, nan or <0 -> most likely one should choose a higer numerical accuracy!!!
    MFP = MFPNanInfOrSmallerZeroCheck(MFP);
    
    return(MFP); 
    
  /*
  if(TotalOrExclusive==1){
    double MFP = tTabulatedTALYSY_Vec[tTabulatedTALYSY_Vec_Component]->GetExclusiveChannelY(InputNucleusId,
											    OutPutNucleusId,
											    LogGamma);
    //check if inf, nan or <0 -> most likely one should choose a higer numerical accuracy!!!
    MFP = MFPNanInfOrSmallerZeroCheck(MFP);
    
    return(MFP);
  }else  if(TotalOrExclusive==2){
    double MFP = tTabulatedTALYSY_Vec[tTabulatedTALYSY_Vec_Component]->GetExclusiveSummedY(InputNucleusId,
											   LogGamma);
    //check if inf, nan or <0 -> most likely one should choose a higer numerical accuracy!!!
    MFP = MFPNanInfOrSmallerZeroCheck(MFP);
    
    return(MFP);
  }
  else{
    throw TCrpErr("Error in TALYSMeanFreePathTabulated::functs(): TotalOrExclusive out of bound."); 
  }
   
   */
}

double TALYSMeanFreePathTabulated::functs(double s, 
					  TVector3D Position,
					  unsigned long int InputNucleusId,
					  unsigned long int OutputId,
					  double NucleusMass,
					  double EnergyOfNucleusGeV,
					  double z,
					  int TotalOrExclusive){
  std::cout<<"Error: functs() shouldn't be called if the tabulated from TALYSMeanFreePathTabulated."<<std::endl;
  exit(-1);
  return(-1.);
}
