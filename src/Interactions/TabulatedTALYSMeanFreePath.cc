/*To calculate the mean free path (MFP) from the TALYS cross section data one has to calculate a double integral which folds the photonbackground (1st integral) and the cross section data (2nd Integral). This can be done using the TALYSMeanFreePath class. Clearly, this is a time consuming procedure which might no be needed for each simulation. In especial, the 2nd integral which mainly averages over the cross section can be tabulated for all simulations. This tables can then be used during the simulation. This tables can be created & loaded by using TabulatedTALYSMeanFreePath.
Aditionally, not all reaction channels are of the same importance. Hence, TabulatedTALYSMeanFreePath can create tables which cointains only those averaged cross sections which are of "importance"; Here that is: Order the reaction channels in ascending order. A reaction channel is "important" if it does not contribute to sum S of all the ordered cross section which cointan e.g. alpha=10% of the total cross section in one energy bin. Note, that the accuracy might be higher than e. g. alpha=10% as the channel is used in the complete energy range if the channel is "important"in only ONE energy bin!!!       
*/

#include "TabulatedTALYSMeanFreePath.h"
#include <algorithm>
#include <TGraph.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TLegend.h>
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
#include "nucleusdb.h"
#include "TALYSMeanFreePathAvrgd.h"

extern std::string TalysDirectory; // from sophiainteractions.cc


TabulatedTALYSMeanFreePath::TabulatedTALYSMeanFreePath(){
  
  fXMin=-1;    
  fXMax=-1;        
  fXNSteps=-1; 
  fXStepsize=-1;
  fthinning_alpha=-1;
  fPDFileType=-1;

  fXUnit="unit free"; 
  fXDescription="#gamma factor of nucleus";
  fYUnit="Mpc^-1";
  fYDescription="inverse mean free path" ;
}

void TabulatedTALYSMeanFreePath::ReadData(){
  TalysDirectory.clear();
  TalysDirectory=fPath.str();
  std::cout<<"Read Files From directory: "<<TalysDirectory<<std::endl;
  
  std::string IntFile = TalysDirectory+"PDExclTabMnFrPthCross.cmt";
  std::ifstream ExclusiveCrossSecData(IntFile.c_str(), std::ios::in);
  if (! ExclusiveCrossSecData.is_open()) 
    throw TCrpErr("Error opening TALYS Cross section data file : " + IntFile
		 +" Note, you have to install them by starting ./GetPDCrossSections.sh"); 
    
  IntFile = TalysDirectory+"PDExclSumTabMnFrPthCross.cmt";
  std::ifstream ExclusiveSumCrossSecData(IntFile.c_str(), std::ios::in);
  if (! ExclusiveSumCrossSecData.is_open()) 
    throw TCrpErr("Error opening TALYS Cross section data file : " + IntFile
		 +" Note, you have to install them by starting ./GetPDCrossSections.sh"); 
  
  IntFile = TalysDirectory+"PDExclTabMnFrPthCrossId.cmt";
  std::ifstream ExclusiveCrossSecDataId(IntFile.c_str(), std::ios::in);
  if (! ExclusiveCrossSecDataId.is_open()) 
    throw TCrpErr("Error opening TALYS Cross section data file : " + IntFile
		 +" Note, you have to install them by starting ./GetPDCrossSections.sh");


  IntFile = TalysDirectory+"PDInitNucleusId.cmt";
  std::ifstream TALYSInitNucleusId(IntFile.c_str(), std::ios::in);
  if (! TALYSInitNucleusId.is_open()) 
    throw TCrpErr("Error opening TALYS Cross section data file : " + IntFile
		 +" Note, you have to install them by starting ./GetPDCrossSections.sh");

  //fill the data from the files into vectors & matrices
  TALYSInitNucleusId 
    >>fPDFileType
    >>fXMin 
    >>fXMax
    >>fXNSteps
    >>fthinning_alpha;
  
   fXStepsize = (fXMax-fXMin) / ( (double) fXNSteps);

   std::cout<<"Minimum gamma factor: "<<fXMin<<std::endl;
   std::cout<<"Maximum gamma factor: "<<fXMax<<std::endl;
   std::cout<<"Number of bins: "<<fXNSteps<<std::endl;
   std::cout<<"thinning factor (100 means no thinning): "<<fthinning_alpha<<std::endl;
   std::cout<<"Read in the mean free path tables for photo-disintegration . . . this may take some time."<<std::endl;
   
  if(fPDFileType!=2){
    std::cout<<"No proper .cmt files were found. This is a big error if you want to run a simulation!"
	     <<" In this case: please, install the cross section package, call make install and restart CRPropa."<<std::endl;
    std::cout<<"It's only a warning if you want to create the .cmt files yourself. Therefore, if you want to create the *.cmt files"
	     <<" yourself please type 1, now."<<std::endl;
    int user;
    std::cin>>user;
    if(user!=1) throw TCrpErr("The files in "+TalysDirectory+" are no proper .cmt files!");
    fPDFileType=1;
  }
  
  double ExclCrossSecVal;
  while (ExclusiveCrossSecData 
	 >>ExclCrossSecVal) {
    fExclusiveYData_Vec.push_back(ExclCrossSecVal) ; 
  }
  
  double ExclSumCrossSecVal;
  while (ExclusiveSumCrossSecData
	 >>ExclSumCrossSecVal) {
    fExclusiveSumYData_Vec.push_back(ExclSumCrossSecVal) ; 
  }
  
  std::vector<unsigned long int> ActualRow(4,0);

   while(ExclusiveCrossSecDataId 
	 >>ActualRow[0] 
	 >>ActualRow[1] 
	 >>ActualRow[2] 
	 >>ActualRow[3]){
     fTALYSExclYId_Matrix.push_back(ActualRow);
   }
   
   std::vector<unsigned long int> ActualRow1x5(5,0);   
   while(TALYSInitNucleusId 
	 >>ActualRow1x5[0] 
	 >>ActualRow1x5[1] 
	 >>ActualRow1x5[2]
	 >>ActualRow1x5[3] 
	 >>ActualRow1x5[4]
	 ){
     fTALYSInitNucleusId_Matrix.push_back(ActualRow1x5);
   }
}

void TabulatedTALYSMeanFreePath::CreateTables(std::vector<PhotonBackground*> PhotonBackgrounds_Vec,
					      std::vector<TabulatedTALYSY*> tTabulatedTALYSY_Vec,
					      PD_TabDataType tPD_TabDataType,
					      TVector3D Position,
					      double epsilon0, 
					      double epsilon_min, 
					      double epsilon_max,
					      std::string folder,
					      int NIntegratorCalls, 
					      double EMinMFPTableVal,
					      double EMaxMFPTableVal,
					      int NBinsMFPTableVal){
 
  if(PhotonBackgrounds_Vec.size()==0 || PhotonBackgrounds_Vec.size()>1){
    throw TCrpErr("TALYSMeanFreePathTabulated::CreateTables(...): No or more than one photonfields were given but can not be handeld.");
  }

  if(tPD_TabDataType != AveragedXSTab 
     &&
     tPD_TabDataType != XSTab) throw TCrpErr("TALYSMeanFreePathTabulated::CreateTables(...): no cross sections or averaged cross sections data given."); 
  
  fFolder=folder;
 
  //open files
  std::string TalysDirectory = folder;

  std::string IntFile = TalysDirectory+"PDExclTabMnFrPthCross.cmt";
  std::ofstream ExclusiveCrossSecData(IntFile.c_str(), std::ios::out);
  if (! ExclusiveCrossSecData.is_open()) 
    throw TCrpErr("Error opening TALYS Cross section data file : " + IntFile); 

  IntFile = TalysDirectory+"PDExclSumTabMnFrPthCross.cmt";
  std::ofstream ExclusiveSumCrossSecData(IntFile.c_str(), std::ios::out);
  if (! ExclusiveSumCrossSecData.is_open()) 
    throw TCrpErr("Error opening TALYS Cross section data file : " + IntFile); 
  
  IntFile = TalysDirectory+"PDExclTabMnFrPthCrossId.cmt";
  std::ofstream ExclusiveCrossSecDataId(IntFile.c_str(), std::ios::out);
  if (! ExclusiveCrossSecDataId.is_open()) 
    throw TCrpErr("Error opening TALYS Cross section data file : " + IntFile);

  IntFile = TalysDirectory+"PDInitNucleusId.cmt";
  std::ofstream TALYSInitNucleusId(IntFile.c_str(), std::ios::out);
  if (! TALYSInitNucleusId.is_open())
    throw TCrpErr("Error opening TALYS Cross section data file : " + IntFile);

  TALYSMeanFreePathAvrgd MeanFreePathAvrgd(epsilon0, epsilon_min, epsilon_max);

  MeanFreePathAvrgd.SetNIntegratorCalls(NIntegratorCalls);

  fXMin=EMinMFPTableVal;
  fXMax=EMaxMFPTableVal;
  fXNSteps=NBinsMFPTableVal;
  fXStepsize = (fXMax-fXMin) / ( (double) fXNSteps);
  
  std::cout<<"CreateTables() for MFP: Set:\nIntegrator Calls\t="<<NIntegratorCalls
	   <<"\nfXMin[MeV]\t="<<fXMin
	   <<"\nfXMax[MeV]\t="<<fXMax
	   <<"\nfXNSteps\t="<<fXNSteps
	   <<"\nfXStepsize[MeV]\t="<<fXStepsize
	   <<std::endl; 
  
  TALYSInitNucleusId<<"2 "<<fXMin<<" "<<fXMax<<" "<<fXNSteps<<" 100"<<std::endl;
  
  std::vector<unsigned long int> ListOfInputNuclei;
  ListOfInputNuclei = tTabulatedTALYSY_Vec[0]->GetListOfInputNuclei();  

  unsigned long int LineCounterOutputNuclei=0;
  unsigned long int LineCounterOutExclChannel=0;
  
  unsigned long int ExclAvrgCrossSecLineCounter=0;
  unsigned long int LineCounterSumExclChannel=0;

  for(int i=0; i < ListOfInputNuclei.size(); i++){
    std::cout<<"process element "<<ListOfInputNuclei[i]<<std::endl;

    TALYSInitNucleusId<<ListOfInputNuclei[i];
    
    TALYSInitNucleusId<<" "<<LineCounterOutExclChannel;  
    
    //////////////////////////////////////////////////////////////////////////////////////////
    //2. Average the exclusive cross sections for the reactions (A,Z)->(#p,#n,#d,#t,#h,#alpha)
    vector<unsigned long int> ListOfAllExclusiveChannels;
    ListOfAllExclusiveChannels=
      tTabulatedTALYSY_Vec[0]
      ->GetListOfAllExclusiveChannels(ListOfInputNuclei[i]);
    
    std::vector<double> sumOfAllExcChannels(fXNSteps, 0.);

    for(int l=0; l < ListOfAllExclusiveChannels.size(); l++){

      //exclude excitation only
      if(ListOfAllExclusiveChannels[l]==0) continue;

      std::cout<<"\t\t -> "<<ListOfAllExclusiveChannels[l]<<std::endl;
      LineCounterOutExclChannel++;
      ExclusiveCrossSecDataId
	<<ListOfAllExclusiveChannels[l]<<" "
	<<ExclAvrgCrossSecLineCounter<<" ";
      
      
      for(int a = 0;              
	  a<fXNSteps;
	  a++){
	
	double logGamma=fXMin+a*(fXMax-fXMin)/(double) fXNSteps;
	
	double FreeMeanPath = MeanFreePathAvrgd.Total_rate_log(PhotonBackgrounds_Vec,
							       tTabulatedTALYSY_Vec,
							       tPD_TabDataType,
							       Position,
							       ListOfInputNuclei[i],
							       ListOfAllExclusiveChannels[l],
							       TNucleusDB::GetInstance()
							       ->GetMass(ListOfInputNuclei[i])*c_squared/1000.,//GeV/c^2
							       pow(10., logGamma)*TNucleusDB::GetInstance()
							       ->GetMass(ListOfInputNuclei[i])*c_squared/1000., //GeV
							       0,                                              //redshift z=0
							       1); //0 means TotalCrossSection
	
	double MFP = FreeMeanPath;
	
	if(1./MFP  > 1.e99)  MFP  = 1./1.e99;
	
	if(ListOfAllExclusiveChannels[l]!=0) sumOfAllExcChannels[a]+=MFP;

	ExclusiveCrossSecData<<MFP<<std::endl;
	ExclAvrgCrossSecLineCounter++;
      }//logGamma
      
      std::vector<int> DMassAndDCharge 
	=  tTabulatedTALYSY_Vec[0]
	->GetDMassAndDChargeExclId(ListOfAllExclusiveChannels[l]);
      ExclusiveCrossSecDataId<<DMassAndDCharge[1]<<" "<<DMassAndDCharge[0]<<std::endl;
      
    }//ListOfExclusiveChannels
    
    //write summed up mean free pathes to file
    for(int b=0; b<sumOfAllExcChannels.size(); b++){
      ExclusiveSumCrossSecData<<sumOfAllExcChannels[b]<<std::endl;
    }

    TALYSInitNucleusId<<" "<<LineCounterOutExclChannel;  
    TALYSInitNucleusId<<" "<<LineCounterSumExclChannel;  
    LineCounterSumExclChannel+=fXNSteps;
    TALYSInitNucleusId<<" "<<LineCounterSumExclChannel;  
    TALYSInitNucleusId<<std::endl;
  }//ListOfInputNuclei
  
  ExclusiveCrossSecData.close();
  ExclusiveSumCrossSecData.close();
  ExclusiveCrossSecDataId.close();
  TALYSInitNucleusId.close();
}


void TabulatedTALYSMeanFreePath::ReduceTables(std::string folder, double accuracy){
  std::cout<<"Start the process to reduce tables with an acuracy of "<<accuracy<<std::endl;

  fthinning_alpha=accuracy;
  
  unsigned long int LineCounterOutputNuclei=0;
  unsigned long int LineCounterOutExclChannel=0;
  unsigned long int LineCounterSumExclChannel=0;
    unsigned long int ExclCrossSecLineCounter=0;
    
  std::string TalysDirectory = folder;

  std::string IntFile = TalysDirectory+"PDTotalTabMnFrPthCross.cmt";

  IntFile = TalysDirectory+"PDExclTabMnFrPthCross.cmt";
  std::ofstream ExclusiveCrossSecData(IntFile.c_str(), std::ios::trunc);
  if (! ExclusiveCrossSecData.is_open()) 
    throw TCrpErr("Error opening TALYS Cross section data file : " + IntFile); 
  
  IntFile = TalysDirectory+"PDExclSumTabMnFrPthCross.cmt";
  std::ofstream ExclusiveSumCrossSecData(IntFile.c_str(), std::ios::out);
  if (! ExclusiveSumCrossSecData.is_open()) 
    throw TCrpErr("Error opening TALYS Cross section data file : " + IntFile);
  
  IntFile = TalysDirectory+"PDExclTabMnFrPthCrossId.cmt";
  std::ofstream ExclusiveCrossSecDataId(IntFile.c_str(), std::ios::trunc);
  if (! ExclusiveCrossSecDataId.is_open()) 
      throw TCrpErr("Error opening TALYS Cross section data file : " + IntFile);

    IntFile = TalysDirectory+"PDInitNucleusId.cmt";
  std::ofstream TALYSInitNucleusId(IntFile.c_str(), std::ios::trunc);
  if (! TALYSInitNucleusId.is_open()) 
    throw TCrpErr("Error opening TALYS Cross section data file : " + IntFile);

  TALYSInitNucleusId<<fPDFileType<<" "
		    <<fXMin<<" "
		    <<fXMax<<" "
		    <<fXNSteps<<" "
		    <<fthinning_alpha<<" "
		    <<std::endl;

  std::vector<unsigned long int> ListOfInputNuclei;
  ListOfInputNuclei = this->GetListOfInputNuclei();  
  
  int AllChannelsCounter=0;
  int FlaggedChannelsCounter=0;
 
  for(int i=0; i < ListOfInputNuclei.size(); i++){
    std::cout<<"process element "<<ListOfInputNuclei[i]<<std::endl;

    vector<unsigned long int> ListOfAllExclusiveChannels;
    ListOfAllExclusiveChannels= this->GetListOfAllExclusiveChannels(ListOfInputNuclei[i]);

    //Calulate the normalization in an energy bin
    int NEnergyBins = this->GetNXBins();
    
    std::vector<  std::vector<unsigned long int>   > * ExclusiveChannelImprtortant;
    ExclusiveChannelImprtortant = new std::vector< std::vector< unsigned long int >  >(ListOfAllExclusiveChannels.size(), std::vector< unsigned long int >(2)); 
    for(int z=0; z<ListOfAllExclusiveChannels.size(); z++){
      (*ExclusiveChannelImprtortant)[z][0]=ListOfAllExclusiveChannels[z];
      (*ExclusiveChannelImprtortant)[z][1]=0;  //0 for false
    }
    
    for(int j=0; j<NEnergyBins; j++){	
      double Normalisation=0.;
      std::vector<  std::vector<double>   > * VectorIdCrossSecMatrix;
      VectorIdCrossSecMatrix = new std::vector< std::vector< double >  >(ListOfAllExclusiveChannels.size(), std::vector< double >(2)); 

      

      for(int l=0; l < ListOfAllExclusiveChannels.size(); l++){  
	double ActVal = this->GetTabExclChannelYValue (ListOfInputNuclei[i], 
						       ListOfAllExclusiveChannels[l], 
						       j);
      
	
	if(  1/ActVal  >= (1.e99 - std::numeric_limits<double>::min())  )   ActVal = 0.;
	
	Normalisation += ActVal;

	(*VectorIdCrossSecMatrix)[l][0] = ActVal;
	(*VectorIdCrossSecMatrix)[l][1] = ListOfAllExclusiveChannels[l];
      }//l: ListOfAllExclusiveChannels
      
      if(Normalisation<=0.){ // => no channel is of importance
	std::cout<<"Normalisation==0 => channels are of no importance"<<std::endl;
	continue;
      }
      
      //sorting of the cross sections
      sort((*VectorIdCrossSecMatrix).begin(), 
	   (*VectorIdCrossSecMatrix).end());
      
      //Check if one of the exclusive channels is of importance in this energy bin
      double ActValue=0.;
      for(int l=0; l < ListOfAllExclusiveChannels.size(); l++){  
	
	double TabVal = this->GetTabExclChannelYValue (ListOfInputNuclei[i], 
						       (long unsigned int) ( (*VectorIdCrossSecMatrix)[l][1]+0.5 ),  
						       j);
	
	if(1./TabVal  >= (1.e99 - std::numeric_limits<double>::min())) {
	  TabVal = 0.;

	  //mark the actual channels as 0 (unimportant)
	  for(int m=0; m<ListOfAllExclusiveChannels.size();m++){
	    if((*ExclusiveChannelImprtortant)[m][0]==ListOfAllExclusiveChannels[l]){
	      (*ExclusiveChannelImprtortant)[m][1]=0;
	      break;
	    }
	    if(m==ListOfAllExclusiveChannels.size()-1){
	      std::cout<<"Warning in TabulatedTALYSMeanFreePath::ReduceTables(...)"
		       <<std::endl;
	      exit(-1);
	    }
	    
	  }//for : m
	  continue;
	}//if
      

	ActValue += TabVal;
	
	if(ActValue*100./Normalisation > 100-accuracy){
	  int ActualExclusiveChannelId = (int) ((*VectorIdCrossSecMatrix)[l][1]+0.5);
	  //make sure to mark the right channel
	  for(int m=0; m<ListOfAllExclusiveChannels.size();m++){
	    if((*ExclusiveChannelImprtortant)[m][0]==ActualExclusiveChannelId){
	      (*ExclusiveChannelImprtortant)[m][1]=1;
	      break;
	    }
	    if(m==ListOfAllExclusiveChannels.size()-1){
	      std::cout<<"Warning in TabulatedTALYSMeanFreePath::ReduceTables(...)"
		       <<std::endl;
	      exit(-1);
	    }
	  }
	}	
      }//l: ListOfAllExclusiveChannels

      delete VectorIdCrossSecMatrix;
    }//j: NEnergyBins
    
    std::cout<<"\n*****RESULTS for "<<ListOfInputNuclei[i]<<std::endl;
    for(int l=0; l < ListOfAllExclusiveChannels.size(); l++){  
      std::cout<<"Is Exclusive channels "<<(*ExclusiveChannelImprtortant)[l][0]
      <<" important? "<<(*ExclusiveChannelImprtortant)[l][1]<<std::endl;
      if(((*ExclusiveChannelImprtortant)[l][1])==1){
	FlaggedChannelsCounter++;
      }
      AllChannelsCounter++;
    }

    FillFiles(folder,
	      ListOfInputNuclei[i],
	      ExclusiveChannelImprtortant,	  
	      &LineCounterOutputNuclei,
	      &LineCounterOutExclChannel,
	      &LineCounterSumExclChannel,
	      &ExclCrossSecLineCounter,
	      &ExclusiveCrossSecData,
	      &ExclusiveCrossSecDataId,
	      &TALYSInitNucleusId,
	      &ExclusiveSumCrossSecData);

    std::cout<<"\n\n"<<std::endl; 
  }//i: ListOfInputNuclei
  
  
  std::cout<<AllChannelsCounter<<" channels tested."
	   <<"\n"<<FlaggedChannelsCounter<<" are of imprtance (accuracy="
	   <<accuracy<<")."<<std::endl;

  ExclusiveCrossSecData.close();
  ExclusiveSumCrossSecData.close();
  ExclusiveCrossSecDataId.close();
  TALYSInitNucleusId.close();
}



void TabulatedTALYSMeanFreePath::FillFiles(std::string folder,
					   unsigned long int  InitNucleusId,
					   std::vector<std::vector <long unsigned int> > * ExclusiveChannelImprtortant,	  
					   unsigned long int* LineCounterOutputNuclei,
					   unsigned long int* LineCounterOutExclChannel,
					   unsigned long int* LineCounterSumExclChannel,
					   unsigned long int* ExclCrossSecLineCounter,
					   ofstream* ExclusiveCrossSecData,
					   ofstream* ExclusiveCrossSecDataId,
					   ofstream* TALYSInitNucleusId,
					   ofstream* ExclusiveSumCrossSecData){
  
  
  (*TALYSInitNucleusId)<<InitNucleusId;
   
  int NEnergyBins = this->GetNXBins();

   (*TALYSInitNucleusId)<<" "<<(*LineCounterOutExclChannel);
  
  
 
 //Exclusive cross sections
  //2. Average the exclusive cross sections for the reactions (A,Z)->(#p,#n,#d,#t,#h,#alpha)
  
  bool GuaranteeOneChannelFlag=false;
  std::vector<double> sumOfAllExcChannels(NEnergyBins, 0.);
  for(int k=0; k < (*ExclusiveChannelImprtortant).size(); k++){
    if((*ExclusiveChannelImprtortant)[k][1]==0) continue;
    GuaranteeOneChannelFlag=true;
    (*LineCounterOutExclChannel)++;
    (*ExclusiveCrossSecDataId)
      <<(*ExclusiveChannelImprtortant)[k][0]<<" "
      <<(*ExclCrossSecLineCounter)<<" ";
    
    for(int l=0; l<NEnergyBins; l++){
      double TabVal = this->GetTabExclChannelYValue (InitNucleusId, 
				   (*ExclusiveChannelImprtortant)[k][0], 
				   l);
      
      double ExclCrossSecVal;
      ExclCrossSecVal= TabVal;

      (*ExclusiveCrossSecData)<<ExclCrossSecVal<<std::endl;
      (*ExclCrossSecLineCounter)++;
      if((*ExclusiveChannelImprtortant)[k][0]!=0.) sumOfAllExcChannels[l]+=ExclCrossSecVal;
    }//l : NEnergyBins  
    
    std::vector<int> DMassAndDCharge 
      =  this->GetDMassAndDChargeExclId((*ExclusiveChannelImprtortant)[k][0]);
    (*ExclusiveCrossSecDataId)<<DMassAndDCharge[1]<<" "<<DMassAndDCharge[0]<<std::endl;
  }// k: ListOfAllExclusiveChannels

   //write summed up mean free pathes to file
  for(int b=0; b<sumOfAllExcChannels.size(); b++){
    (*ExclusiveSumCrossSecData)<<sumOfAllExcChannels[b]<<std::endl;
  }

  if(GuaranteeOneChannelFlag==false){
    std::cerr<<"WARNING: TabulatedTALYSMeanFreePath::FillFiles() for Reduced Tables."
	     <<" No channel has survived choose a less small alpha."
	     <<std::endl;
  }

  (*TALYSInitNucleusId)<<" "<<(*LineCounterOutExclChannel);
  (*TALYSInitNucleusId)<<" "<<(*LineCounterSumExclChannel);  
  (*LineCounterSumExclChannel)+=fXNSteps;
  (*TALYSInitNucleusId)<<" "<<(*LineCounterSumExclChannel)<<std::endl;

}

