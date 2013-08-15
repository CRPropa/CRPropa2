/**
 @file    TabulatedTALYSY.cc
 @author  Nils Nierstenhoefer, nierstenhofer@physik.uni-wuppertal.de
 @brief   Implementation of the TabulatedTALYSY class. See the .h file
*/

#include "TabulatedTALYSY.h"
#include <algorithm>
#include <sstream>
#include <iostream>
#include <fstream>

TabulatedTALYSY::TabulatedTALYSY() 
  : fXMin(-1) 
  , fXMax(-1)        
  , fXNSteps(-1) 
  , fXStepsize(-1.)
 {}

 TabulatedTALYSY::~TabulatedTALYSY() {
  
}

double TabulatedTALYSY::GetTabExclChannelYValue(unsigned long int InputNucleusId,
						unsigned long int ExclusiveChannelId,
						int EnergyBin){
  
  //firstly identify the input nucleus via NucleusId = Chargenumber*1000+MassNumber 
  unsigned long int ExclusiveYDataId_Beg=0;
  unsigned long int ExclusiveYDataId_End=0;
  
  std::vector< std::vector< unsigned long int >  >::iterator the_iterator;
  
  bool found=false;
  the_iterator=fTALYSInitNucleusId_Matrix.begin();
  while( the_iterator != fTALYSInitNucleusId_Matrix.end() && !found ) {
    if((*the_iterator)[0]==InputNucleusId){
      found=true;
      break;
    }
    the_iterator++;
  }
  if(found==false){
    return(YValueNotAvailable());
  }
  

  ExclusiveYDataId_Beg = (*the_iterator)[1];
  ExclusiveYDataId_End = (*the_iterator)[2];
  
  found=false;
  unsigned long int ExclusiveYData_Beg=0;
  for(long unsigned int i=ExclusiveYDataId_Beg; i<ExclusiveYDataId_End; i++){
    if(fTALYSExclYId_Matrix[i][0]==ExclusiveChannelId) {
      ExclusiveYData_Beg=fTALYSExclYId_Matrix[i][1];
      found=true;
      break;    
    }
  }
  
  
  if(found==false){
     std::cerr<<" WARNING  : There is no exclusive channel crossection for "<<InputNucleusId<<" -> "
	      <<ExclusiveChannelId<<std::endl;
     return(YValueNotAvailable());
  }
  return(fExclusiveYData_Vec[ExclusiveYData_Beg+EnergyBin]);
}




double TabulatedTALYSY::GetTabExclSummedYValue(unsigned long int InputNucleusId,
					       int EnergyBin){  
  //firstly identify the input nucleus via NucleusId = Chargenumber*1000+MassNumber 
  unsigned long int ExclusiveSumYDataId_Beg=0;
  unsigned long int ExclusiveSumYDataId_End=0;
   
  std::vector< std::vector< unsigned long int >  >::iterator the_iterator;
  
  bool found=false;
  the_iterator=fTALYSInitNucleusId_Matrix.begin();
  while( the_iterator != fTALYSInitNucleusId_Matrix.end() ) {
    if((*the_iterator)[0]==InputNucleusId){
      found=true;
      break;
    }
    the_iterator++;
  }
  if(found==false){
    //std::cerr<<"Warning: There is no summed, exclusive crossection for InputNucleusId="
    //<<InputNucleusId<<std::endl;
    return(YValueNotAvailable());
  }
  
  ExclusiveSumYDataId_Beg = (*the_iterator)[3];
  ExclusiveSumYDataId_End = (*the_iterator)[4];
  
  //sanity check for testing purpose only : TODO: comment out. 
  if(ExclusiveSumYDataId_Beg+EnergyBin>ExclusiveSumYDataId_End){
     throw TCrpErr("Error: in TabulatedTALYSY::GetTabExclSummedYValue(...): ExclusiveSumYDataId_Beg+EnergyBin>ExclusiveSumYDataId_End");
  }

  return(fExclusiveSumYData_Vec[ExclusiveSumYDataId_Beg+EnergyBin]);
}





int TabulatedTALYSY::GetEnergyBinIdBelow(double Energy){
  if(Energy<fXMin ||  Energy>fXMax){
    //std::cerr<<"ERROR 4 in TabulatedTALYSY::GetEnergyBinIdBelow: Energy="<<Energy
    //     <<" out of range: "
    //     <<"("<<fXMin<<","<<fXMax<<"] MeV"<<std::endl;
    
    //exit(-1);
  }
  return( (int) ((Energy-fXMin)/fXStepsize));
}


 int TabulatedTALYSY::GetEnergyBinIdAbove(double Energy){
  return(GetEnergyBinIdBelow(Energy)+1);
}

 double TabulatedTALYSY::GetEnergyForEnergyBinId(unsigned long int EnergyBin){
  return(fXMin+EnergyBin*1.*fXStepsize);
}


double TabulatedTALYSY::GetExclusiveChannelY(unsigned long int InputNucleusId,
					     int ExclusiveChannelId,
					     double Energy){
  if(Energy<=GetEnergyForEnergyBinId(0)){
    return(ExclYValueBelowFXMin(InputNucleusId, ExclusiveChannelId));
  }
   
   if(Energy >= GetEnergyForEnergyBinId(GetNXBins()-1)){
     return(ExclYValueAboveFXMax(InputNucleusId, ExclusiveChannelId));
   }


   //Get EnergyBinId of the bin below/above.
   int EnBinIdBel = GetEnergyBinIdBelow(Energy);
   int EnBinIdAbv = GetEnergyBinIdAbove(Energy);
      
   //DeltaEn
   double EnergyBelow=GetEnergyForEnergyBinId(EnBinIdBel);
   double EnergyAbove=GetEnergyForEnergyBinId(EnBinIdAbv);
   double DeltaEn    =EnergyAbove-EnergyBelow;
   
   //DeltaExclY
   double ExclYBel = GetTabExclChannelYValue(InputNucleusId,
					     ExclusiveChannelId,
					     EnBinIdBel);
   double ExclYAbv = GetTabExclChannelYValue(InputNucleusId,
					     ExclusiveChannelId,
					     EnBinIdAbv);
   double DeltaExclY      =ExclYAbv-ExclYBel;

   return(ExclYBel+(DeltaExclY/DeltaEn)*(Energy-EnergyBelow));
 }





double TabulatedTALYSY:: GetExclusiveSummedY(unsigned long int InputNucleusId,
					     double Energy){
  
   if(Energy<=GetEnergyForEnergyBinId(0)){
     return(ExclSumYValueBelowFXMin(InputNucleusId));
   }
   
   if(Energy >= GetEnergyForEnergyBinId(GetNXBins()-1)){
     return(ExclSumYValueAboveFXMax(InputNucleusId));
   }


   //Get EnergyBinId of the bin below/above.
   int EnBinIdBel = GetEnergyBinIdBelow(Energy);
   int EnBinIdAbv = GetEnergyBinIdAbove(Energy);
   
   //DeltaEn
   double EnergyBelow=GetEnergyForEnergyBinId(EnBinIdBel);
   double EnergyAbove=GetEnergyForEnergyBinId(EnBinIdAbv);
   double DeltaEn    =EnergyAbove-EnergyBelow;
   
   //DeltaExclY
   double ExclSumYBel        = GetTabExclSummedYValue(InputNucleusId,
						      EnBinIdBel);
   double ExclSumYAbv        = GetTabExclSummedYValue(InputNucleusId,
						      EnBinIdAbv);
   double DeltaExclSumY      = ExclSumYAbv-ExclSumYBel;
   
  return(ExclSumYBel+(DeltaExclSumY/DeltaEn)*(Energy-EnergyBelow));
}


 void TabulatedTALYSY::RootPlotExclY(unsigned long int InputNucleusId,
			      unsigned long int ExclusiveChannelId){

  std::vector<double> energy;
  std::vector<double> cross_section;
  std::vector<double> discrete_values;
  int counter=0;
  
  for(double e=fXMin; e<=fXMax;e+=(fXMax-fXMin)/50000.){
    double value = GetExclusiveChannelY(InputNucleusId,
					ExclusiveChannelId,
					e);
    energy.push_back(e);
    cross_section.push_back(value);
    counter++;
    
    int BinBelow = GetEnergyBinIdBelow(e);
    discrete_values.push_back(
			      GetTabExclChannelYValue(InputNucleusId,
						      ExclusiveChannelId,
						      BinBelow)
			      );
    
  }  
  
  TGraph* ContinousAproxY=new TGraph(counter, 
				     &(energy[0]), 
				     &(cross_section[0]));
  
  TGraph* DiscreteAproxY=new TGraph(counter, 
				    &(energy[0]), 
				    &(discrete_values[0]));
  
  
  std::ostringstream title;
  title<<"NucelusTransition+InputNucleusId:"
       <<InputNucleusId
       <<"->"
       <<ExclusiveChannelId;
  
  TCanvas* Can = new TCanvas(title.str().c_str(),title.str().c_str(),1);
  Can->SetTitle(title.str().c_str());
  
  ContinousAproxY->Draw("Alp");
  ContinousAproxY->SetTitle(title.str().c_str());
  ContinousAproxY->SetMarkerColor(2);
  ContinousAproxY->SetLineColor(2);
  ContinousAproxY->GetXaxis()->SetTitle((fXDescription+" ["+fXUnit+"]").c_str());
  std::ostringstream reaction;
  reaction<<"Total:"<<InputNucleusId
	  <<"->"
	  <<ExclusiveChannelId;
  ContinousAproxY->GetYaxis()->SetTitle((fYDescription+"("+reaction.str()+") ["+fYUnit+"]").c_str());
  
  DiscreteAproxY->Draw("SAMEpl");
   DiscreteAproxY->GetXaxis()->SetTitle((fXDescription+" ["+fXUnit+"]").c_str());
   DiscreteAproxY->GetYaxis()->SetTitle((fYDescription+"("+reaction.str()+") ["+fYUnit+"]").c_str());

   TLegend* legend = new TLegend(.75,.80,.95,.95);
   legend->AddEntry(ContinousAproxY,"Continous aprox. of cross section");
   legend->AddEntry(DiscreteAproxY,"discrete cross section values");
   legend->Draw();


   std::ostringstream filename;
   filename<<"Exclusive_"<<InputNucleusId
	   <<"_xs"
	   <<ExclusiveChannelId;
   Can->SaveAs((filename.str()+".root").c_str());
   delete Can;
   delete ContinousAproxY;
   delete DiscreteAproxY;
   delete legend;
}

void TabulatedTALYSY::ASCIIOutExclY(unsigned long int InputNucleusId,
				    unsigned long int ExclusiveChannelId){
  std::ostringstream name;
  name<<"ExclusiveY"<<InputNucleusId<<"_xs"<<ExclusiveChannelId<<".txt";
  
  ofstream fp_out(name.str().c_str(), ios::trunc);
  
  fp_out<<"\n\t\t**Exclusive "<<fYDescription<<" for :"
	<<InputNucleusId<<"-> xs"<<ExclusiveChannelId<<std::endl;
  fp_out<<"\t\t"<<fYDescription<<" ["<<fXUnit<<"]\t\tY ["<<fYUnit<<"]"<<std::endl;
  for(int i=0; i<=fXNSteps-1; i++){
    fp_out<<"\t\t"<<GetEnergyForEnergyBinId(i)
	  <<"\t\t"<<GetTabExclChannelYValue(InputNucleusId, ExclusiveChannelId, i)
	  <<std::endl;
  }
  fp_out.close();
}

std::vector<unsigned long int>  
TabulatedTALYSY::GetListOfInputNuclei(){

  std::vector< std::vector< unsigned long int >  >::iterator the_iterator;
  
  std::vector<long unsigned int> ListOfInputNucleiVector;
  //find the input nucleus
  the_iterator=fTALYSInitNucleusId_Matrix.begin();
  while( the_iterator != fTALYSInitNucleusId_Matrix.end() ) {
    ListOfInputNucleiVector.push_back((*the_iterator)[0]);
    the_iterator++;
  }
  return(ListOfInputNucleiVector);
}

 std::vector<int> TabulatedTALYSY::GetDMassAndDChargeExclId (int ExclusiveChannelId){

  //#neutrons
   int NNeutron         = (int) (ExclusiveChannelId/1.e5);
   ExclusiveChannelId -= (long unsigned int) (NNeutron*1.e5);

  //#protons
   int NProton         = (int) (ExclusiveChannelId/1.e4);
   ExclusiveChannelId-= (long unsigned int) (NProton*1.e4);

  //#deuterium
   int NDeuterium      = (int) (ExclusiveChannelId/1.e3);
   ExclusiveChannelId-= (long unsigned int) (NDeuterium*1.e3);

  //#tritium
   int NTritium        = (int) (ExclusiveChannelId/1.e2);
   ExclusiveChannelId-= (long unsigned int) (NTritium*1.e2);

  //#Helium
   int NHelium         = (int) (ExclusiveChannelId/1.e1);
   ExclusiveChannelId-= (long unsigned int) (NHelium*1.e1);

   int NAlpha          = (int) (ExclusiveChannelId);
   
   int DeltaMass=NNeutron+NProton+2*NDeuterium+3*(NTritium+NHelium)+4*NAlpha;
   int DeltaCharge=NProton+NDeuterium+NTritium+2*NHelium+2*NAlpha;
   
   std::vector<int> DMassAndDChargeVec;
   
   DMassAndDChargeVec.push_back(DeltaMass);
   DMassAndDChargeVec.push_back(DeltaCharge);

  return(DMassAndDChargeVec);

}  


std::vector<unsigned long int> 
TabulatedTALYSY::GetMassAndChargeNuclId(unsigned long int NucleusId){

  unsigned long int ChargeNumber = (unsigned long int) (NucleusId/1000.);
  unsigned long int MassNumber   = (unsigned long int) (NucleusId-ChargeNumber*1000.);

  std::vector<unsigned long int> GetDMassAndDChargeNuclIdVec; 
  GetDMassAndDChargeNuclIdVec.push_back(MassNumber);
  GetDMassAndDChargeNuclIdVec.push_back(ChargeNumber);
  
  return(GetDMassAndDChargeNuclIdVec);
}


std::vector<int> 
TabulatedTALYSY::GetDMassAndDChargeNuclId(unsigned long int NucleusId1,
					  unsigned long int NucleusId2){
  
  
  std::vector<unsigned long int> GetDMassAndDChargeNuclIdVec_1;
  GetDMassAndDChargeNuclIdVec_1=this->GetMassAndChargeNuclId(NucleusId1);
  
  std::vector<unsigned long int> GetDMassAndDChargeNuclIdVec_2;
  GetDMassAndDChargeNuclIdVec_2=this->GetMassAndChargeNuclId(NucleusId2);

  std::vector<int> DMassAndDChargeNuclIdVec;
  DMassAndDChargeNuclIdVec.push_back(GetDMassAndDChargeNuclIdVec_1[0]
				     -
				     GetDMassAndDChargeNuclIdVec_2[0]);

  DMassAndDChargeNuclIdVec.push_back(GetDMassAndDChargeNuclIdVec_1[1]
				     -
				     GetDMassAndDChargeNuclIdVec_2[1]);

  return(DMassAndDChargeNuclIdVec);
}




vector<unsigned long int> 
TabulatedTALYSY::GetListOfAllExclusiveChannels(unsigned long int InputNucleusId){
  //firstly identify the input nucleus via NucleusId = Chargenumber*1000+MassNumber 
   unsigned long int ExclusiveYDataId_Beg=0;
   unsigned long int ExclusiveYDataId_End=0;
   
   vector<unsigned long int> OutputExclVec;

   std::vector< std::vector< unsigned long int >  >::iterator the_iterator;
   
   bool found=false;
   the_iterator=fTALYSInitNucleusId_Matrix.begin();
   while( the_iterator != fTALYSInitNucleusId_Matrix.end() ) {
     if((*the_iterator)[0]==InputNucleusId){
       found=true;
       break;
     }
     the_iterator++;
   }

   if(found==false){
     return(OutputExclVec);
   }
   
   ExclusiveYDataId_Beg = (*the_iterator)[1];
   ExclusiveYDataId_End = (*the_iterator)[2];
   
   found=false;
   for(long unsigned int i=ExclusiveYDataId_Beg; i<ExclusiveYDataId_End; i++){
     OutputExclVec.push_back(fTALYSExclYId_Matrix[i][0]);
     found=true;
   }
   if(found==false){
     //std::cerr<<"ERROR 7: There is no Exclusive cross section for Id="
     //	       <<InputNucleusId
     //	       <<std::endl;
      //exit(-1);
    }
   return(OutputExclVec);
}










  double test=0.;

TGraph* TabulatedTALYSY::RootPlotYVsX(unsigned long int InputNucleusId,
				      long unsigned int OutId,
				      double fXMin,
				      double fXMax,
				      double fXNSteps,
				      double redshift,
				      int TotalSumOrExclusive){ 

  std::vector<double> EnVec;
  std::vector<double> CalculatedMeanFreePathVec;
  double MeanFreePathSum     = 0.;
  double MeanFreePathSumZZero= 0.;

  ofstream fp_out;

  std::vector<double> XVec;
  std::vector<double> YVec;

  for(double e=fXMin; e<=fXMax; 
      e += (fXMax-fXMin)/((double) (fXNSteps))){
    XVec.push_back(e);
    double TheYVal;
  
    if(TotalSumOrExclusive==1) TheYVal = this->GetExclusiveChannelY(InputNucleusId, OutId, e);
    
    if(TotalSumOrExclusive==2) {
      TheYVal = this->GetExclusiveSummedY(InputNucleusId, e);
    }
    YVec.push_back(TheYVal);
  }

  if(TotalSumOrExclusive==1){
    test+=GetExclusiveChannelY(InputNucleusId, OutId, fXMax);
  }
  TGraph* TheGraph = new TGraph(XVec.size(), &(XVec[0]), &(YVec[0]));
  ostringstream Title;
  Title<<InputNucleusId<<"_"<<OutId;
  TheGraph->SetName(Title.str().c_str()); 
  TheGraph->SetTitle(Title.str().c_str()); 
  TheGraph->GetXaxis()->SetTitle((fXDescription+" ["+fXUnit+"]").c_str());
  TheGraph->GetYaxis()->SetTitle((fYDescription+" ["+fYUnit+"]").c_str());
  return( TheGraph );
}


TGraph* TabulatedTALYSY::RootPlotExclSumYVsX(unsigned long int InputNucleusId,
					     double fXMin,
					     double fXMax,
					     double fXNSteps,
					     double redshift){
 std::vector<double> EnVec;
 std::vector<double> CalculatedMeanFreePathVec;

 std::vector<double> XVec;
 std::vector<double> YVec;
 int counter=0;

 for(double e=fXMin; e<=fXMax; 
     e += (fXMax-fXMin)/((double) (fXNSteps))){
   
   XVec.push_back(e);
   YVec.push_back(0);

   
   std::vector<unsigned long int> OutPutNucleiVec = this->GetListOfAllExclusiveChannels(InputNucleusId);

   for(int o=0; o < OutPutNucleiVec.size(); o++){

     if(OutPutNucleiVec[o]==0) continue;     
     std::vector<int> MassChargeVec = this->GetDMassAndDChargeExclId(OutPutNucleiVec[o]);   	
     
     double NucleusMass =  TNucleusDB::GetInstance()->GetMass(InputNucleusId)*c_squared/1000.;
     double TheYVal = this->GetExclusiveChannelY(InputNucleusId,
						 OutPutNucleiVec[o],
						 e);
     YVec[counter]+=TheYVal; 
   }
   counter++;
 }
  TGraph* TheGraph = new TGraph(XVec.size(), &(XVec[0]), &(YVec[0]));
  ostringstream Title;
  Title<<InputNucleusId;
  TheGraph->SetName(Title.str().c_str()); 
  TheGraph->SetTitle(Title.str().c_str()); 
  TheGraph->GetXaxis()->SetTitle("#gamma energy/MeV");
  TheGraph->GetYaxis()->SetTitle("#sigma/mb");
  return( TheGraph );
}


void TabulatedTALYSY::RootFileWithAllReactions(double fXMin,
					       double fXMax,
					       double fXNSteps,
					       double redshift){
  
  std::cout<<"void TSophiaInteractions::RootFileWithAllReactions(double fXMin, double fXMax, double fXNSteps, double redshift) called"
	   <<std::endl;

  TFile* theRootFileExcl            = new TFile("AllInOneFileExcl.root","RECREATE");
  TFile* theRootFileExclSummed      = new TFile("AllInOneFileExclSummed.root","RECREATE");
  TFile* theRootFileExclSummedTable = new TFile("AllInOneFileExclSummedTable.root","RECREATE");
  
  std::vector< unsigned long int > InPutNucleiVec = this->GetListOfInputNuclei();
  
  //exclusive cross sections
  for(int i=0; i < InPutNucleiVec.size(); i++){
    std::cout<<"\n\n*****exclsive (single)\nInputNucleusId="<<InPutNucleiVec[i]<<std::endl;
    
    std::vector<unsigned long int> OutPutNucleiVec =
      this->GetListOfAllExclusiveChannels(InPutNucleiVec[i]);
    for(int o=0; o < OutPutNucleiVec.size(); o++){
      theRootFileExcl->cd();
      
      ostringstream Title;
      Title<<InPutNucleiVec[i]<<"_"<<OutPutNucleiVec[o]<<std::endl;
      
      TCanvas* TheCan = new TCanvas(("Can"+Title.str()).c_str(),("Can"+Title.str()).c_str(),1); 
      TGraph* TheGraph = RootPlotYVsX(InPutNucleiVec[i],
				      OutPutNucleiVec[o],
				      fXMin,
				      fXMax,
				      fXNSteps,
				      redshift,
				      1);
       
       TheGraph->Draw("A*");
       TheGraph->Write();
       delete TheGraph;
       delete TheCan;
    }
  }
  theRootFileExcl->Close();
  delete theRootFileExcl;
  //exclusive cross sections summed up to total MFP
   for(int i=0; i < InPutNucleiVec.size(); i++){
     std::cout<<"\n\n*****exclsive (sum)\nInputNucleusId="<<InPutNucleiVec[i]<<std::endl;
     theRootFileExclSummed->cd();
     TGraph* TheGraph = RootPlotExclSumYVsX(InPutNucleiVec[i],				      
					    fXMin,
					    fXMax,
					    fXNSteps,
					    redshift);
     
     TheGraph->Draw("A*");
     TheGraph->Write();
     TheGraph->GetXaxis()->SetTitle("#gamma energy/MeV");
     TheGraph->GetYaxis()->SetTitle("#sigma/mb");
     delete TheGraph; 
   }
   theRootFileExclSummed->Close();
   delete theRootFileExclSummed;
}

