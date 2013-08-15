/*To calculate the mean free path (MFP) from the TALYS cross section data one has to calculate a double integral which folds the photonbackground (1st integral) and the cross section data (2nd Integral). This can be done using the TALYSMeanFreePath class. Clearly, this is a time consuming procedure which might no be needed for each simulation. In especial, the 2nd integral which mainly averages over the cross section can be tabulated for all simulations. This tables can then be used during the simulation. This tables can be created & loaded by using TabulatedTALYSAveragedCrossSection.
Aditionally, not all reaction channels are of the same importance. Hence, TabulatedTALYSAveragedCrossSection can create tables which cointains only those averaged cross sections which are of "importance"; Here that is: Order the reaction channels in ascending order. A reaction channel is "important" if it does not contribute to sum S of all the ordered cross section which cointan e.g. alpha=10% of the total cross section in one energy bin. Note, that the accuracy might be higher than e. g. alpha=10% as the channel is used in the complete energy range if the channel is "important"in only ONE energy bin!!!       
*/

#include "TabulatedTALYSAveragedCrossSection.h"

TabulatedTALYSAveragedCrossSection* TabulatedTALYSAveragedCrossSection::finstance;

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
#include "TALYSMeanFreePathAccurate.h"

TabulatedTALYSAveragedCrossSection::TabulatedTALYSAveragedCrossSection(){
  
  fPD_TabDataType=AveragedXSTab;
  fXUnit="MeV"; 
  fXDescription="upper border inner integral";
  fYUnit="mbarn*GeV*GeV";
  fYDescription="avarged cross section (inner integral)" ;

  //open files
  std::string TalysDirectory = DEFAULT_TALYSAveraged_DIR;
  
  std::cout<<"Read Files From directory: "<<TalysDirectory<<std::endl;

  std::string IntFile = "";

  IntFile = TalysDirectory+"PDExclAvrgdCross.cax";
  std::ifstream ExclusiveCrossSecData(IntFile.c_str(), std::ios::in);
  if (! ExclusiveCrossSecData.is_open()) 
  throw TCrpErr("Error opening TALYS Cross section data file : " + IntFile
		 +" Note, you have to install them by starting ./GetPDCrossSections.sh"); 
  
  IntFile = TalysDirectory+"PDExclSumAvrgdCross.cax";
  std::ifstream ExclusiveSumCrossSecData(IntFile.c_str(), std::ios::in);
  if (! ExclusiveSumCrossSecData.is_open()) 
    throw TCrpErr("Error opening TALYS Cross section data file : " + IntFile
		 +" Note, you have to install them by starting ./GetPDCrossSections.sh"); 
    
  IntFile = TalysDirectory+"PDExclAvrgdCrossId.cax";
  std::ifstream ExclusiveCrossSecDataId(IntFile.c_str(), std::ios::in);
  if (! ExclusiveCrossSecDataId.is_open()) 
    throw TCrpErr("Error opening TALYS Cross section data file : " + IntFile
		 +" Note, you have to install them by starting ./GetPDCrossSections.sh");
  
    IntFile = TalysDirectory+"PDInitNucleusId.cax";
  std::ifstream TALYSInitNucleusId(IntFile.c_str(), std::ios::in);
  if (! TALYSInitNucleusId.is_open()) 
    throw TCrpErr("Error opening TALYS Cross section data file : " + IntFile
		 +" Note, you have to install them by starting ./GetPDCrossSections.sh");

  TALYSInitNucleusId 
    >>fPDFileType
    >>fXMin 
    >>fXMax
    >>fXNSteps
    >>fthinning_alpha;

  fXStepsize = (fXMax-fXMin) / ( (double) fXNSteps);
  
  if(fPDFileType!=1){
    std::cout<<"No proper .cax files were found. This is a big error if you want to run a simulation!"
	     <<" In this case: please, install the cross section package, call make install and restart CRPropa."<<std::endl;
    std::cout<<"It's only a warning if you want to create the .cax files yourself. Therefore, if you want to create the *.cax files"
	     <<" yourself please type 1, now."<<std::endl;
    int user;
    std::cin>>user;
    if(user!=1) throw TCrpErr("The files in "+TalysDirectory+" are no proper .cxs files!");
    fPDFileType=1;
  }

  std::cout<<"Minimum gamma energy [MeV]: "<<fXMin<<std::endl;
  std::cout<<"Maximum gamma energy [MeV]: "<<fXMax<<std::endl;
  std::cout<<"Number of bins: "<<fXNSteps<<std::endl;
  std::cout<<"Thinning-factor (100 means no thinnig): "<<fthinning_alpha<<std::endl;
  std::cout<<"Read in the averaged cross section data for photo-disintegration . . . this may take some time."<<std::endl;
  //fill the data from the files into vectors & matrices

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
  //matrices (created using the std:vector class)
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

TabulatedTALYSAveragedCrossSection::~TabulatedTALYSAveragedCrossSection() {
  
}

void TabulatedTALYSAveragedCrossSection::CreateTables(std::string folder){
  std::cout<<"void TabulatedTALYSAveragedCrossSection::CreateTables(. . . )"<<std::endl;

  fFolder=folder;
 
  //open files
  std::string TalysDirectory = folder;

  std::cout<<"Write Files to:"<<TalysDirectory<<std::endl;

  std::string IntFile = "";
  
  IntFile = TalysDirectory+"PDExclAvrgdCross.cax";
  std::ofstream ExclusiveCrossSecData(IntFile.c_str(), std::ios::out);
  if (! ExclusiveCrossSecData.is_open()) 
  throw TCrpErr("Error opening TALYS Cross section data file : " + IntFile); 
 
   IntFile = TalysDirectory+"PDExclSumAvrgdCross.cax";
  std::ofstream ExclusiveSumCrossSecData(IntFile.c_str(), std::ios::out);
  if (! ExclusiveSumCrossSecData.is_open()) 
    throw TCrpErr("Error opening TALYS Cross section data file : " + IntFile); 
 
  IntFile = TalysDirectory+"PDExclAvrgdCrossId.cax";
  std::ofstream ExclusiveCrossSecDataId(IntFile.c_str(), std::ios::out);
  if (! ExclusiveCrossSecDataId.is_open()) 
  throw TCrpErr("Error opening TALYS Cross section data file : " + IntFile);

  IntFile = TalysDirectory+"PDInitNucleusId.cax";
  std::ofstream TALYSInitNucleusId(IntFile.c_str(), std::ios::out);
  if (! TALYSInitNucleusId.is_open())
    throw TCrpErr("Error opening TALYS Cross section data file : " + IntFile);
  

  std::cout<<"1. Give the minimum energy for table (e.g. 0.001 MeV):"<<std::endl;
  std::cin>>fXMin;
  std::cout<<"2. Give the maxmimum energy for table (e.g. 250 MeV):"<<std::endl;
  std::cin>>fXMax;
  std::cout<<"3. Give number of bins (integer) for the tables (e.g. N=500):"<<std::endl;
  std::cin>>fXNSteps;

  fXStepsize = (fXMax-fXMin) / ( (double) fXNSteps);

  TALYSInitNucleusId<<"1 "<<fXMin<<" "<<fXMax<<" "<<fXNSteps<<" 100"<<std::endl;

  std::vector<unsigned long int> ListOfInputNuclei;
  ListOfInputNuclei = TabulatedTALYSCrossSection::GetInstance()->GetListOfInputNuclei();  
  
  TALYSMeanFreePathAccurate MeanFreePathAccurate;
  
  unsigned long int LineCounterOutputNuclei=0;
  unsigned long int LineCounterOutExclChannel=0;
  
  unsigned long int ExclAvrgCrossSecLineCounter=0;
  unsigned long int LineCounterSumExclChannel=0;

  for(int i=0; i < ListOfInputNuclei.size(); i++){
    std::cout<<"process element "<<ListOfInputNuclei[i]<<std::endl;
    
    TALYSInitNucleusId<<ListOfInputNuclei[i]<<" ";
    TALYSMeanFreePathAccurate TALYSMFPAcc;
    TALYSInitNucleusId<<" "<<LineCounterOutExclChannel;  
    
    //////////////////////////////////////////////////////////////////////////////////////////
    //2. Average the exclusive cross sections for the reactions (A,Z)->(#p,#n,#d,#t,#h,#alpha)
    vector<unsigned long int> ListOfAllExclusiveChannels;
    ListOfAllExclusiveChannels=
      TabulatedTALYSCrossSection::GetInstance()
      ->GetListOfAllExclusiveChannels(ListOfInputNuclei[i]);
    
    std::vector<double> sumOfAllExcChannels(fXNSteps, 0.);

    for(int l=0; l < ListOfAllExclusiveChannels.size(); l++){
      
      //exclude excitation only
      if(ListOfAllExclusiveChannels[l]==0) continue;

      LineCounterOutExclChannel++;
      
      ExclusiveCrossSecDataId
	<<ListOfAllExclusiveChannels[l]<<" "
	<<ExclAvrgCrossSecLineCounter<<" ";

      for(int a = 0;              
	  a<fXNSteps;
	  a++){
	
	double IntegralUpperBorder=fXMin+a*(fXMax-fXMin)/(double) fXNSteps;
	
	//There is no position dependence here . . . position is needed because of coding procedures
	TVector3D Position;
	Position.setX(0);
	Position.setY(0);
	Position.setZ(0);
	double MFP = TALYSMFPAcc.functs_int(IntegralUpperBorder/1000.,
					    Position,
					    ListOfInputNuclei[i],
					    ListOfAllExclusiveChannels[l],
					    TNucleusDB::GetInstance()
					    ->GetMass(ListOfInputNuclei[i])*c_squared/1000.,
					    0.,
					    0.,
					    1);  //1 means exclusive cross section
	ExclusiveCrossSecData<<MFP<<std::endl;
	
	if(ListOfAllExclusiveChannels[l]!=0) sumOfAllExcChannels[a]+=MFP;

	ExclAvrgCrossSecLineCounter++;
      }//IntegralUpperBorder
      
      std::vector<int> DMassAndDCharge 
	=  TabulatedTALYSCrossSection::GetInstance()
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


void TabulatedTALYSAveragedCrossSection::ReduceTables(std::string folder, double accuracy){

  std::cout<<"void TabulatedTALYSAveragedCrossSection::ReduceTables(. . . )"<<std::endl;
  
  fthinning_alpha=accuracy;
  
  unsigned long int LineCounterOutputNuclei=0;
  unsigned long int LineCounterOutExclChannel=0;

  unsigned long int ExclCrossSecLineCounter=0;
  unsigned long int LineCounterSumExclChannel=0;
    
  std::string TalysDirectory = folder;

  std::string IntFile = "";

  IntFile = TalysDirectory+"PDExclAvrgdCross.cax";
  std::ofstream ExclusiveCrossSecData(IntFile.c_str(), std::ios::trunc);
  if (! ExclusiveCrossSecData.is_open()) 
    throw TCrpErr("Error opening TALYS Cross section data file : " + IntFile); 
    
  IntFile = TalysDirectory+"PDExclSumAvrgdCross.cax";
  std::ofstream ExclusiveSumCrossSecData(IntFile.c_str(), std::ios::trunc);
  if (! ExclusiveSumCrossSecData.is_open()) 
    throw TCrpErr("Error opening TALYS Cross section data file : " + IntFile); 
  
  IntFile = TalysDirectory+"PDExclAvrgdCrossId.cax";
  std::ofstream ExclusiveCrossSecDataId(IntFile.c_str(), std::ios::trunc);
  if (! ExclusiveCrossSecDataId.is_open()) 
    throw TCrpErr("Error opening TALYS Cross section data file : " + IntFile);

    IntFile = TalysDirectory+"PDInitNucleusId.cax";
  std::ofstream TALYSInitNucleusId(IntFile.c_str(), std::ios::trunc);
  if (! TALYSInitNucleusId.is_open()) 
    throw TCrpErr("Error opening TALYS Cross section data file : " + IntFile);

  TALYSInitNucleusId<<"1 "<<fXMin<<" "<<fXMax<<" "<<fXNSteps<<" "<<fthinning_alpha<<std::endl;

  std::vector<unsigned long int> ListOfInputNuclei;
  ListOfInputNuclei = TabulatedTALYSAveragedCrossSection::GetInstance()->GetListOfInputNuclei();  
  
  int AllChannelsCounter=0;
  int FlaggedChannelsCounter=0;
 
  for(int i=0; i < ListOfInputNuclei.size(); i++){
    std::cout<<"process element "<<ListOfInputNuclei[i]<<std::endl;

    vector<unsigned long int> ListOfAllExclusiveChannels;
    ListOfAllExclusiveChannels= TabulatedTALYSAveragedCrossSection::GetInstance()
      ->GetListOfAllExclusiveChannels(ListOfInputNuclei[i]);

    //Calulate the normalization an energy bin
    int NEnergyBins = TabulatedTALYSAveragedCrossSection::GetInstance()->GetNXBins();
    
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
	double ActVal = TabulatedTALYSAveragedCrossSection::GetInstance()
	  ->GetTabExclChannelYValue (ListOfInputNuclei[i], 
				     ListOfAllExclusiveChannels[l], 
				     j);
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
      
	//Check if one of the exclusive channels of importance in this energy bin
      double ActValue=0.;
      for(int l=0; l < ListOfAllExclusiveChannels.size(); l++){  
        ActValue += TabulatedTALYSAveragedCrossSection::GetInstance()
	  ->GetTabExclChannelYValue (ListOfInputNuclei[i], 
				     (int) ( (*VectorIdCrossSecMatrix)[l][1]+0.5 ), 
				     j);
	
	if(ActValue*100./Normalisation > 100-accuracy){
	  int ActualExclusiveChannelId = (int) ((*VectorIdCrossSecMatrix)[l][1]+0.5);
	  //make sure to mark the right channel
	  for(int m=0; m<ListOfAllExclusiveChannels.size();m++){
	    if((*ExclusiveChannelImprtortant)[m][0]==ActualExclusiveChannelId){
	      (*ExclusiveChannelImprtortant)[m][1]=1;
	      break;
	    }
	    if(m==ListOfAllExclusiveChannels.size()-1){
	      std::cout<<"Error in TabulatedTALYSAveragedCrossSection::ReduceTables(...)"
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



void TabulatedTALYSAveragedCrossSection::FillFiles(
			   std::string folder,
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
  
  std::cout<<"Start: void TabulatedTALYSAveragedCrossSection::FillFiles( . . . )"<<std::endl;
  
  (*TALYSInitNucleusId)<<InitNucleusId;

  int NEnergyBins = TabulatedTALYSAveragedCrossSection::GetInstance()->GetNXBins();
  
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
      double ExclCrossSecVal = TabulatedTALYSAveragedCrossSection::GetInstance()
	->GetTabExclChannelYValue (InitNucleusId, 
				   (*ExclusiveChannelImprtortant)[k][0], 
				   l);
      
      (*ExclusiveCrossSecData)<<ExclCrossSecVal<<std::endl;
      (*ExclCrossSecLineCounter)++;
      if((*ExclusiveChannelImprtortant)[k][0]!=0.) sumOfAllExcChannels[l]+=ExclCrossSecVal;
    }//l : NEnergyBins  
    
    std::vector<int> DMassAndDCharge 
      =  TabulatedTALYSAveragedCrossSection::GetInstance()
      ->GetDMassAndDChargeExclId((*ExclusiveChannelImprtortant)[k][0]);
    (*ExclusiveCrossSecDataId)<<DMassAndDCharge[1]<<" "<<DMassAndDCharge[0]<<std::endl;
  }// k: ListOfAllExclusiveChannels

   //write summed up mean free pathes to file
  for(int b=0; b<sumOfAllExcChannels.size(); b++){
    (*ExclusiveSumCrossSecData)<<sumOfAllExcChannels[b]<<std::endl;
  }

  
  
  (*TALYSInitNucleusId)<<" "<<(*LineCounterOutExclChannel);
  (*TALYSInitNucleusId)<<" "<<(*LineCounterSumExclChannel);  
  (*LineCounterSumExclChannel)+=fXNSteps;
  (*TALYSInitNucleusId)<<" "<<(*LineCounterSumExclChannel)<<std::endl;
}


