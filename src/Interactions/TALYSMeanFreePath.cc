#include "TALYSMeanFreePath.h"
#include "TALYSMeanFreePathTabulated.h"
#include "TabulatedTALYSMeanFreePath.h"
#include "PhotonBackground.h"
#include <TCanvas.h>
#include <TGraph.h>
#include <ostream>
#include <sstream>
#include <limits>

TALYSMeanFreePath::TALYSMeanFreePath(double epsilon_zero, double epsilon_min, double epsilon_max)
: _fEpsilon0(epsilon_zero), _fEpsilon_min(epsilon_min), _fEpsilon_max(epsilon_max){
}

TALYSMeanFreePath::TALYSMeanFreePath(){
}

TALYSMeanFreePath::~TALYSMeanFreePath(){}

double PHOTD(double t);

//RETURNS total interaction rate with CMB in Mpc^-1
double TALYSMeanFreePath::Total_rate(std::vector<PhotonBackground*> PhotonBkgrd,
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
  
  if(PhotonBkgrd.size()==0) throw TCrpErr("Total_rate(...): No Photonfields were given!"); 
  if(tTabulatedTALYSY_Vec.size()>1) throw TCrpErr("Total_rate(...): tTabulatedTALYSY_Vec.size()>1!"); 
  if(tPD_TabDataType != AveragedXSTab 
     && 
     tPD_TabDataType != XSTab) throw TCrpErr("Total_rate(...): Got no cross sections or averaged cross sections data."); 

  double  xmp = NucleusMass;          
  double Pp = sqrt(EnergyOfNucleusGeV*EnergyOfNucleusGeV-xmp*xmp);  
  
  int numberOfGaussIntegratorCalls=_fNIntegrationCalls; 
  double result=0.;
  
  ClassMemPtr Ptr;
  Ptr=&TALYSMeanFreePath::functs_int_times_PhtnBckGrnd;

  double gamma = EnergyOfNucleusGeV/NucleusMass;

  //loop over photon backgrounds
  for(int m=0; m<PhotonBkgrd.size(); m++){

    fPhotonBkgrd = PhotonBkgrd[m];
    //Todo: This is a temoprary solution as CMB/IRB mean free path is treated by one class. Instead, one sould derive 
    //TALYSMeanFreePath*_IRB and TALYSMeanFreePath*_CMB and set _fepsilon0 via the corresponding contrcutor!
    _fEpsilon0   = fPhotonBkgrd->GetEpsilonZero(); 

    //Guenters suggested integration range
    double UpperBorder = 9.99999920000006370e-01;
    double lowerBorder = 4.03225806451612898e-22; //ICRCPoster
  
    for(int i=0; i<=numberOfGaussIntegratorCalls-1; i++){
      
      //Borders if equally binned in the space of the substitued (E=kT(1-t)/t) variable
      double lower=lowerBorder+(UpperBorder-lowerBorder)*i/(double) numberOfGaussIntegratorCalls;
      double upper=lowerBorder+(UpperBorder-lowerBorder)*(i+1)/(double) numberOfGaussIntegratorCalls;
      
      result+=gauss(Ptr, 
		    Position,
		    lower,
		    upper,
		    InputNucleusId,
		    OutPutNucleusId,
		    NucleusMass,
		    EnergyOfNucleusGeV,
		    0.,        //set redshift to zero according to formula (13) in CRPropa paper 
		    TotalOrExclusive
		    );
    }// i : number of gauss integrator calls  
  }// m : photon fields under consideration
  
  result=1.302e16*3.086e22/2./(EnergyOfNucleusGeV/NucleusMass)/(EnergyOfNucleusGeV/NucleusMass)*result;
  
  //check if inf, nan or <0 -> most likely one should choose a higer numerical accuracy!!!
  result = MFPNanInfOrSmallerZeroCheck(result);
  return(result);
  
}

void TALYSMeanFreePath::SetNIntegratorCalls(int NIntegrationCalls){
  _fNIntegrationCalls=NIntegrationCalls;
}

 /*This member calculates the MFP on the fly with subst. t=ln(E/epsilon0) in the outer integral*/
double TALYSMeanFreePath::Total_rate_log(std::vector<PhotonBackground*> PhotonBkgrd,
					 std::vector<TabulatedTALYSY*> tTabulatedTALYSY_Vec,
					 PD_TabDataType tPD_TabDataType,
					 TVector3D Position,
					 unsigned long int InputNucleusId,
					 unsigned long int OutPutNucleusId,
					 double NucleusMass,
					 double EnergyOfNucleusGeV,
					 double z,
					 int TotalOrExclusive){
  
  if(PhotonBkgrd.size()==0) throw TCrpErr("Total_rate(...): No Photonfields were given!"); 
  if(tTabulatedTALYSY_Vec.size()>1) throw TCrpErr("Total_rate(...): tTabulatedTALYSY_Vec.size()>1!"); 
  if(tPD_TabDataType != AveragedXSTab 
     && 
     tPD_TabDataType != XSTab) throw TCrpErr("Total_rate(...): Got no cross sections or averaged cross sections data."); 


  double  xmp = NucleusMass;          
  double Pp = sqrt(EnergyOfNucleusGeV*EnergyOfNucleusGeV-xmp*xmp);  
  
  int numberOfGaussIntegratorCalls=_fNIntegrationCalls; //ICRCPoster
  
  double result=0.;
  
  ClassMemPtr Ptr;
  Ptr=&TALYSMeanFreePath::functs_int_times_PhtnBckGrnd_log;
  
  double gamma = EnergyOfNucleusGeV/NucleusMass;

  //loop over photon backgrounds
  for(int m=0; m<PhotonBkgrd.size(); m++){
    fPhotonBkgrd = PhotonBkgrd[m];
    _fEpsilon0   = fPhotonBkgrd->GetEpsilonZero();   
    
    //Guenter's suggestions (new ln)
    double lowerBorder = log(_fEpsilon_min / _fEpsilon0);
    double UpperBorder = log(_fEpsilon_max / _fEpsilon0);
    
    for(int i=0; i<=numberOfGaussIntegratorCalls-1; i++){
      
      double lower=lowerBorder+(UpperBorder-lowerBorder)*i/(double) numberOfGaussIntegratorCalls;
      double upper=lowerBorder+(UpperBorder-lowerBorder)*(i+1)/(double) numberOfGaussIntegratorCalls;
      
      result+=gauss(Ptr,
		    Position,
		    lower,
		    upper,
		    InputNucleusId,
		    OutPutNucleusId,
		    NucleusMass,
		    EnergyOfNucleusGeV,
		    0.,        //set redshift to zero according to formula (13) in CRPropa paper 
		    TotalOrExclusive
		    );

    } // i : number of gauss integrator calls  
  } // m : photon fields under consideration
  
  result=1.302e16*3.086e22/2./(EnergyOfNucleusGeV/NucleusMass)/(EnergyOfNucleusGeV/NucleusMass)*result;
  
  //check if inf, nan or <0 -> most likely one should choose a higer numerical accuracy!!!
  result = MFPNanInfOrSmallerZeroCheck(result);
  return(result);
}


double TALYSMeanFreePath::MFPNanInfOrSmallerZeroCheck(double MFP){

  if(isnan(MFP)){
    std::cout<<"Warning: MFP is nan! Consider to increase the numerical accuracy"<<std::endl;
    std::cout<<"Set to std::numeric_limits<double>::min()"<<std::endl;
    return(std::numeric_limits<double>::min());
  }
  if(isinf(MFP)){
    std::cout<<"Warning: MFP is inf! Consider to increase the numerical accuracy"<<std::endl;
    std::cout<<"Set to std::numeric_limits<double>::min()"<<std::endl;
    return(std::numeric_limits<double>::min());
  }
  if(MFP<=0.){
    //std::cout<<"Warning: MFP is <=0! Consider to increase the numerical accuracy"<<std::endl;
    //std::cout<<"Set to 1.e-99"<<std::endl;
    return(std::numeric_limits<double>::min());
  }
  

  return(MFP);
}

double TALYSMeanFreePath::functs_int_times_PhtnBckGrnd(double t,
						       TVector3D Position,
						       unsigned long int InputNucleusId,
						       unsigned long int OutPutNucleusId,
						       double NucleusMass,
						       double EnergyOfNucleusGeV,
						       double z,
						       int TotalOrExclusive){

  double epsilon0 = _fEpsilon0; 
    
  double gamma = EnergyOfNucleusGeV/NucleusMass;
  double UpperBorder2ndInt = 2*gamma*epsilon0*(1-t)/t; 

  double PBG_density = fPhotonBkgrd->GetPhotonDensity(Position.getX()/Mpc, 
						      Position.getY()/Mpc, 
						      Position.getZ()/Mpc,  
						      z, 
						      epsilon0*(1-t)/t);

  double AveragedCrossSection=functs_int(UpperBorder2ndInt,
					 Position,
					 InputNucleusId,
					 OutPutNucleusId,
					 NucleusMass,
					 EnergyOfNucleusGeV,
					 z,
					 TotalOrExclusive);
  return(PBG_density/(1-t)/(1-t)/epsilon0*AveragedCrossSection);         
}



 double  TALYSMeanFreePath::functs_int_times_PhtnBckGrnd_log(double t,
							     TVector3D Position,
							     unsigned long int InputNucleusId,
							     unsigned long int OutPutNucleusId,
							     double NucleusMass,
							     double EnergyOfNucleusGeV,
							     double z,
							     int TotalOrExclusive){
   double epsilon0 = _fEpsilon0; 
   
   

   double gamma = EnergyOfNucleusGeV/NucleusMass;
   double UpperBorder2ndInt = 2*gamma*exp(t)*_fEpsilon0; 
   double PBG_density = fPhotonBkgrd->GetPhotonDensity(Position.getX()/Mpc, 
						       Position.getY()/Mpc, 
						       Position.getZ()/Mpc, 
						       z, 
						       exp(t)*_fEpsilon0);

   double AveragedCrossSection=functs_int(UpperBorder2ndInt,
					  Position,
					  InputNucleusId,
					  OutPutNucleusId,
					  NucleusMass,
					  EnergyOfNucleusGeV,
					  z,
					  TotalOrExclusive);

   return(PBG_density / _fEpsilon0 / exp(t) * AveragedCrossSection);
 }
	  


//double testUpperBorder=0;
//integrand for Total_rate_cmb / inner of the double integral
double TALYSMeanFreePath::functs_int(double UpperBorder,
				     TVector3D Position,
				     unsigned long int InputNucleusId,
				     unsigned long int OutPutNucleusId,
				     double NucleusMass,
				     double EnergyOfNucleusGeV,
				     double z,
				     int TotalOrExclusive){
 
  double resultSum=0.;
  double resultGauss=0.;
  
  if (UpperBorder < 0.) {
    resultSum=0.;
    resultGauss=0.;
  }
  else{
    ClassMemPtr Ptr;
    Ptr=&TALYSMeanFreePath::functs;
    
    //to have a constant ammount of sampling point in cross section data for different upper borders in the inner integral.
    double stepsize = fabs(250-1.e-3)/((double) 500);
    
    int NBins = (int) ( ((UpperBorder*1000-1.e-3)/stepsize)/8. );
    if(NBins==0) NBins=1;
    
    int NumberOfGaussIntegratorCalls=NBins;
    for(int j=0; j<NumberOfGaussIntegratorCalls ; j++){
      resultGauss+=gauss(Ptr,
			 Position,
			 j*UpperBorder/((double) NumberOfGaussIntegratorCalls),
			 (j+1.)*UpperBorder/((double) NumberOfGaussIntegratorCalls),
			 InputNucleusId,
			 OutPutNucleusId,
			 NucleusMass,
			 EnergyOfNucleusGeV,
			 z,
			 TotalOrExclusive);
    }
  }   
  return(resultGauss);          //Z ?
}

//Photonbackground:
//= 0 (CMB only)
//= 1 (IRO only)
//= 2 (CMB+IRO)
//= 3 (Tabulated (assuming constant photonfield))
TGraph* TALYSMeanFreePath::RootPlotMFPVsX(std::vector<PhotonBackground*> PhotonBkgrd,
					  std::vector<TabulatedTALYSY*> tTabulatedTALYSY_Vec,
					  TVector3D Position,
					  unsigned long int InputNucleusId,
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
    double NucleusMass =  TNucleusDB::GetInstance()->GetMass(InputNucleusId)*c_squared/1000.;
    
    XVec.push_back(e);
    
    double FreeMeanPath=0;

    for(int k=0; k<tTabulatedTALYSY_Vec.size(); k++){
      
      FreeMeanPath += pow(1+redshift,3)*this->Total_rate_log(PhotonBkgrd,
							     tTabulatedTALYSY_Vec,
							     tTabulatedTALYSY_Vec[k]->GetPD_TabDataType(),
							     Position,
							     InputNucleusId,
							     OutId,
							     NucleusMass,
							     pow(10,e)*NucleusMass*(1+redshift),
							     0.,       //s. eq 13 from CRPropa paper
							     TotalSumOrExclusive); //0 means TotalCrossSection
      
    }
    FreeMeanPath=1/FreeMeanPath;
    YVec.push_back(FreeMeanPath);
  }
  
  TGraph* TheGraph = new TGraph(XVec.size(), &(XVec[0]), &(YVec[0]));
  ostringstream Title;
  Title<<InputNucleusId<<"_"<<OutId<<"_"<<TotalSumOrExclusive;
  TheGraph->SetName(Title.str().c_str()); 
  TheGraph->SetTitle(Title.str().c_str()); 
  return( TheGraph );
}

//Photonbackground:
//= 0 (CMB only)
//= 1 (IRO only)
//= 2 (CMB+IRO)
//= 3 (Tabulated (assuming constant photonfield))
TGraph* TALYSMeanFreePath::RootPlotExclSumMFPVsX(std::vector<PhotonBackground*> PhotonBkgrd,
						 std::vector<TabulatedTALYSY*> tTabulatedTALYSY_Vec,
						 TVector3D Position,
						 unsigned long int InputNucleusId,
						 double fXMin,
						 double fXMax,
						 double fXNSteps,
						 double redshift){
 std::vector<double> EnVec;
 std::vector<double> CalculatedMeanFreePathVec;
 double MeanFreePathSum = 0.;

 std::vector<double> XVec;
 std::vector<double> YVec;
 int counter=0;
 
 //create common list of all exclusive channels from IRB and CMB
 std::vector<unsigned long int> AllExclusiveChannels =
   this->GetCombinedExclusiveChannelList(tTabulatedTALYSY_Vec,
					 InputNucleusId);
 
 for(double e=fXMin; e<=fXMax; 
     e += (fXMax-fXMin)/((double) (fXNSteps))){
   
   XVec.push_back(e);
   YVec.push_back(0);  

   for(int o=0; o < AllExclusiveChannels.size(); o++){
     
     if(AllExclusiveChannels[o]==0) continue;
     double NucleusMass =  TNucleusDB::GetInstance()->GetMass(InputNucleusId)*c_squared/1000.;

     for(int k=0; k<tTabulatedTALYSY_Vec.size(); k++){
       double FreeMeanPath=0;

      FreeMeanPath = pow(1+redshift,3)*this->Total_rate_log(PhotonBkgrd,
							     tTabulatedTALYSY_Vec,
							     tTabulatedTALYSY_Vec[k]->GetPD_TabDataType(),
							     Position,
							     InputNucleusId,
							     AllExclusiveChannels[o],
							     NucleusMass,
							     pow(10,e)*NucleusMass*(1+redshift),
							     0.,       //s. eq 13 from CRPropa paper
							     1); //0 means TotalCrossSection
       
       YVec[counter]+=FreeMeanPath;
     }
   }
   
   YVec[counter]=1./YVec[counter];
   counter++;
 }
  TGraph* TheGraph = new TGraph(XVec.size(), &(XVec[0]), &(YVec[0]));
  ostringstream Title;
  Title<<InputNucleusId<<"_3";
  TheGraph->SetName(Title.str().c_str()); 
  TheGraph->SetTitle(Title.str().c_str()); 
  return( TheGraph );
}






void TALYSMeanFreePath::RootFileWithAllMFPs(std::vector<PhotonBackground*> PhotonBkgrd,
					    std::vector<TabulatedTALYSY*> tTabulatedTALYSY_Vec,
					    TVector3D Position,
					    double fXMin,
					    double fXMax,
					    double fXNSteps,
					    double redshift){
  
  std::cout<<"void TALYSMeanFreePath::RootFileWithAllMFPs(double fXMin, double fXMax, double fXNSteps, double redshift) called"
	   <<std::endl;
  
  TFile* theRootFileExcl = new TFile("AllMFPsInOneFileExcl.root","RECREATE");
  TFile* theRootFileExclSummed = new TFile("AllMFPsInOneFileExclSummed.root","RECREATE");
  TFile* theRootFileExclSummedTable = new TFile("AllMFPsInOneFileExclSummedTable.root","RECREATE");

  std::vector< unsigned long int > InPutNucleiVec = tTabulatedTALYSY_Vec[0]->GetListOfInputNuclei();

   //exclusive cross sections
  for(int i=0; i < InPutNucleiVec.size(); i++){
    std::cout<<"\n\n*****exclsive (single)\nInputNucleusId="<<InPutNucleiVec[i]<<std::endl;
    std::vector<unsigned long int> OutPutNucleiVec =
      this->GetCombinedExclusiveChannelList(tTabulatedTALYSY_Vec,
					    InPutNucleiVec[i]);

    for(int o=0; o < OutPutNucleiVec.size(); o++){
      theRootFileExcl->cd();
      
      if(OutPutNucleiVec[o]==0) continue;

      ostringstream Title;
      Title<<InPutNucleiVec[i]<<"_"<<OutPutNucleiVec[o]<<"_1"<<std::endl;
      
       TCanvas* TheCan = new TCanvas(("Can"+Title.str()).c_str(),("Can"+Title.str()).c_str(),1); 
       TGraph* TheGraph = RootPlotMFPVsX(PhotonBkgrd,
					 tTabulatedTALYSY_Vec,
					 Position,
					 InPutNucleiVec[i],
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
    TGraph* TheGraph = RootPlotExclSumMFPVsX(PhotonBkgrd,
					     tTabulatedTALYSY_Vec,
					     Position,
					     InPutNucleiVec[i],				      
					     fXMin,
					     fXMax,
					     fXNSteps,
					     redshift);
       
       TheGraph->Draw("A*");
       TheGraph->Write();
       delete TheGraph;
   }
   theRootFileExclSummed->Close();
   delete theRootFileExclSummed;

   
   //exclusive cross sections summed up to total MFP
   for(int i=0; i < InPutNucleiVec.size(); i++){
     std::cout<<"\n\n*****exclsive (sum from table)\nInputNucleusId="<<InPutNucleiVec[i]<<std::endl;
     theRootFileExclSummedTable->cd();
     TGraph* TheGraph = RootPlotMFPVsX(PhotonBkgrd,
				       tTabulatedTALYSY_Vec,
				       Position,
				       InPutNucleiVec[i],
				       0,
				       fXMin,
				       fXMax,
				       fXNSteps,
				       redshift,
				       2);
     
     TheGraph->Draw("A*");
     TheGraph->Write();
     delete TheGraph;
   }
   theRootFileExclSummedTable->Close();
   delete theRootFileExclSummedTable;
}

std::vector<unsigned long int> TALYSMeanFreePath::GetCombinedExclusiveChannelList(std::vector<TabulatedTALYSY*> tTabulatedTALYSY_Vec,
									   unsigned long int NucleusId){
  //the one component case is simple
  if(tTabulatedTALYSY_Vec.size()==1) return tTabulatedTALYSY_Vec[0]->GetListOfAllExclusiveChannels(NucleusId);
  
  //more then one component case
  std::vector<unsigned long int> buffer;
  for(int k=0; k<tTabulatedTALYSY_Vec.size(); k++){
    std::vector<unsigned long int> next = tTabulatedTALYSY_Vec[k]->GetListOfAllExclusiveChannels(NucleusId);;
    std::vector<unsigned long int> old  = buffer;
    buffer.clear();
    buffer = GetUniqueVecFrom2Vec(next, old);
  }
  
  return buffer;
}


std::vector<unsigned long int> TALYSMeanFreePath::GetUniqueVecFrom2Vec(std::vector<unsigned long int> vec1,
										std::vector<unsigned long int> vec2){
  std::vector<unsigned long int> buffer(vec1.size()+vec2.size());
    if(vec1.size()==0) return vec2;
    if(vec2.size()==0) return vec1;
    if(vec1.size()==0 && vec2.size()==0) return buffer;

  sort(vec1.begin(), vec1.end());
  sort(vec2.begin(), vec2.end());

  merge(vec1.begin(), vec1.end(),
	vec2.begin(), vec2.end(),
	buffer.begin());

  sort(buffer.begin(), buffer.end());
  
  std::vector<unsigned long int>::iterator it;
  it = unique(buffer.begin(), buffer.end());
  
  buffer.resize( it - buffer.begin() );

 return(buffer);
}
