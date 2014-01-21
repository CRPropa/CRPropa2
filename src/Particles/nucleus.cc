/**
 @file    nucleus.cc
 @author  Eric Armengaud, armengau@in2p3.fr
 @brief   Implementation of the TNucleus class. See the .h file
*/

#include "nucleus.h"
#include "sophia.h"
#include "sibyll.h"
#include "CRPropa.h"

#include "TALYSMeanFreePathAccurate.h"
#include "TALYSMeanFreePathAvrgd.h"
#include "TALYSMeanFreePathTabulated.h"
#include "TabulatedTALYSAveragedCrossSection.h"
#include "TIRBzEvolutionModel.h"

#include <cstdio>
#include <ctime>  
#include <limits>
#include <sys/time.h>


double TestDistanceSum_Excl=0.;


extern long pprod, bbarprod, nprod; 
extern double interactt, deflectt, photodt, pionprodt,  pairpt, decayt;

TNucleus::TNucleus(int aMassNumber, int aChargeNumber)  {
    SetType(PARTICLE_NUCLEUS) ;
  _fMassNumber = aMassNumber ;
  _fChargeNumber = aChargeNumber ; 
  _fMass = this->Mass() ;
  _fCharge = eplus *_fChargeNumber ;
  _fTimeStep = 1 * Mpc * inv_c_light ;
  _fNextTimeStep = _fTimeStep ;
  _fpUniverse = NULL;
  _fInitType=1000 * aChargeNumber + aMassNumber;
  throw  TCrpErr("This constructor should not be used!" ); 
  }

TNucleus::TNucleus(TUniverse* aUniv, 
		   TVector3D aPosition, 
		   TVector3D aMomentum, 
		   TVector3D aInitPosition, 
		   TVector3D aInitMomentum, 
		   double aTime, 
		   int aMassNumber,
		   int aChargeNumber,
		   TList1DPhotons* apList1DPhotons, 
		   int aInitType
		   ) {
  // Constructor useful when several nuclei are generated
  
  //Test
  //if(aChargeNumber==0 && aMassNumber==1) std::cout<<"neutron creted"<<std::endl;
  
  SetType(PARTICLE_NUCLEUS);
  _fpList1DPhotons = apList1DPhotons ;
  _fMassNumber     = aMassNumber ;
  _fChargeNumber   = aChargeNumber ; 
  _fMass           = this->Mass() ;
  _fCharge         = eplus *_fChargeNumber ;
  _fTime           = aTime ;
  _fPosition       = aPosition ;
  _fMomentum       = aMomentum ;
  _fInitPosition   = aInitPosition ;
  _fInitMomentum   = aInitMomentum ;
  _fEnergy         = c_light * _fMomentum.mag() ;
  _fpUniverse      = aUniv;
  _fInitType       = aInitType;
    da = _fpUniverse->DistanceArray();
    
  if(aMassNumber < aChargeNumber)  throw  TCrpErr("Secondary constructor: aMassNumber < aChargeNumber!" );  
  if((_fPosition.x() < _fpUniverse->Xmin()) || (_fPosition.x() > _fpUniverse->Xmax())){
    std::cout << "_fPosition.x()" << _fPosition.x()/Mpc <<std::endl;
    std::cout << "_fpUniverse->Xmax()" << _fpUniverse->Xmax()/Mpc <<std::endl;
    std::cout << " _fpUniverse->Xmin()" <<  _fpUniverse->Xmin()/Mpc <<std::endl;
    throw TCrpErr("Particle outside box X @ beginning of constructor" );
  }
  if((_fPosition.y() < _fpUniverse->Ymin()) || (_fPosition.y() > _fpUniverse->Ymax()))
    throw TCrpErr("Particle outside box Y @ beginning of contructor" );
  if((_fPosition.z() < _fpUniverse->Zmin()) || (_fPosition.z() > _fpUniverse->Zmax()))
    throw TCrpErr("Particle outside box Z @ beginning of contructor" );

  // Timesteps and distance to observers
  _fDistObs = 1 * Gpc ;
  if (aUniv->Observers()->Type() == OBSERVER_SMALLSPHERE ) {
    TObservers* lObs = aUniv->Observers() ;
    _fWasInsideSpheres.assign(aUniv->Observers()->Nb(),false) ;
    /*std::cout<<"Create & initialize _fWasInsideSpheres for small sphere case with (constructor 2)"
	     <<_fWasInsideSpheres.size()<<" elements."
	     <<std::endl;*/
    for (int j=0; j<lObs->Nb(); j++) {
      double lDistObs = (_fPosition - lObs->Positions(j)).mag() - lObs->Radius() ;
      //check if particle was created inside the observer; in this case set _fWasInsideSpheres=true!
      if (lDistObs > 0) {
	_fDistObs = min(_fDistObs,lDistObs) ;
      }
      else{
	//std::cout<<"Particle created in sphere. Set _fWasInsideSpheres=true (constructor 2)"<<std::endl;
	_fWasInsideSpheres[j]=true;	
      }
    }
  }
  if (aUniv->Observers()->Type() == OBSERVER_LARGESPHERE ){
    //We need to know if secondary particle is inside or outside of the initial sphere.
    for (int j=0; j<_fpUniverse->Observers()->Nb(); j++) {
      double lDistObs = (_fPosition-_fInitPosition).mag() - _fpUniverse->Observers()->Radii(j) ;
      if(lDistObs < 0){      
	_fIsInsideSpheres.assign(aUniv->Observers()->Nb(),1) ;
      } else {
	_fIsInsideSpheres.assign(aUniv->Observers()->Nb(),0) ;
      }
    }
  }
  //_fTimeStep = 1.e-6 * min(aUniv->InteractionData()->InteractionTimeStep(),_fDistObs* inv_c_light) ;
  _fTimeStep = aUniv->IntegratorMinTimeStep() ;
  _fNextTimeStep = _fTimeStep ;
  
  // Initial redshift
  if (aUniv->Type() != UNIVERSE_ENV1D  || !_fpUniverse->InteractionData()->RedshiftFlag()) {
    _fRedshift = 0 ; // will be used for rate computations by default
  } else{
    this->ComputeRedshift() ;
    _fInitRedshift=_fRedshift;
    //std::cout<<"redshift="<<_fInitRedshift<<std::endl;
    //    std::cout<<"A="<<_fMassNumber<<std::endl;
    //    std::cout<<"Z="<<_fChargeNumber<<std::endl;
  }
  
  if((_fPosition.x() < _fpUniverse->Xmin()) || (_fPosition.x() > _fpUniverse->Xmax()))
    throw TCrpErr("Particle outside box X @ end of constructor" );
  if((_fPosition.y() < _fpUniverse->Ymin()) || (_fPosition.y() > _fpUniverse->Ymax()))
    throw TCrpErr("Particle outside box Y @ end of constructor" );
  if((_fPosition.z() < _fpUniverse->Zmin()) || (_fPosition.z() > _fpUniverse->Zmax()))
    throw TCrpErr("Particle outside box Z @ end of constructor" );  

}

//constructor used for particle injection
TNucleus::TNucleus(TUniverse* aUniv, 
		   TList1DPhotons* apList1DPhotons) {
  SetType(PARTICLE_NUCLEUS);

  
  TestDistanceSum_Excl=0.;

  TSources* sources = aUniv->Sources();
  _sDetected=0;
  _fTime           = 0 ;
  _fpList1DPhotons = apList1DPhotons ;
  _fpUniverse      = aUniv;
    da = _fpUniverse->DistanceArray();

  //_fcounter=0;

  //std::cout << "New trajectory" << std::endl; 
  int Amax=1;
  if(sources->RigidityFlag()){
    for(int i=0; i < sources->GetInitialMassNumber().size(); i++){
      Amax=(Amax < sources->GetInitialMassNumber().at(i))? sources->GetInitialMassNumber().at(i) : Amax;
    }
  }// else Amax=1;


#ifdef DEBUG_OUTPUT
  static int flag=1;
  if(flag){flag=0; std::cout << " Amax " << Amax << " sources->RigidityFlag() " << sources->RigidityFlag()  << std::endl;}
#endif

  do {


#ifdef DEBUG_OUTPUT
	std::cout << "!!!!!Start of propagate!!!" << std::endl;
#endif

    // 0 ) Mass and Charge
    int i=0;
    double proba=RandFlat::shoot();

    while(sources->GetNucleiAbundance().at(i) < proba) i++; 
    _fMassNumber     = sources->GetInitialMassNumber()[i];
    _fChargeNumber   = sources->GetInitialChargeNumber()[i]; 
    _fMass           = this->Mass() ;
    _fCharge         = eplus *_fChargeNumber ; 
    _fInitType=1000 * _fChargeNumber + _fMassNumber;
    
    // 1 ) Initial position and energy
    //  TSources *sources = _fpUniverse->Sources() ;    
    
    if(_fMassNumber < _fChargeNumber)  throw  TCrpErr("Injection constructor: _fMassNumber < _fChargeNumber!" );  
    if (sources->Type() == SOURCE_DISCRETE) {
      // Pick up a source from the list
      int i = RandFlat::shootInt(sources->Nb()) ;
      _fPosition = sources->Positions(i) ;
      
      if ( sources->SpectrumFlag() == 1) { // monochromatic spectrum
	if( sources->RigidityFlag()) {
	  _fEnergy = sources->EcutList(i) * double(_fChargeNumber);
	} else _fEnergy = sources->EcutList(i);
      } else if (sources->SpectrumFlag() == 2) { // spectral index
	double toto = RandFlat::shoot() ;
	double alpha = sources->AlphaList(i) ;
	double lEmax = sources->EcutList(i) ;
	if (sources->RigidityFlag()==1) lEmax*=double(_fChargeNumber);
	double lEmin = sources->EminList(i) ;

// 	if ( alpha == 1 ) {
// 	  _fEnergy = pow(sources->Emin()/double(Amax),toto)*pow(lEmax,1-toto) ;
// 	} else {
// 	  _fEnergy = pow(pow(sources->Emin()/double(Amax),1-alpha)+
// 			 toto*(pow(lEmax,1-alpha) -pow(sources->Emin()/double(Amax),1-alpha)),1/(1-alpha)) ;
// 	}
	_fEnergy=RandPowerLaw(-alpha, lEmin,lEmax,toto);

      }
    } else if (sources->Type() == SOURCE_CONTINUOUS) {
      _fPosition = sources->GetPosition() ;
      
      if ( sources->SpectrumFlag() == 1) { // monochr
	if( sources->RigidityFlag()) {
	  _fEnergy = sources->Ecut() * double(_fChargeNumber);
	} else _fEnergy = sources->Ecut();
      } else if (sources->SpectrumFlag() == 2) { // spectral index
	double toto = RandFlat::shoot() ;
	double alpha = sources->Alpha() ;
	double lEmax = sources->Ecut() ;
	if (sources->RigidityFlag()==1) lEmax*=double(_fChargeNumber);
// 	if ( alpha == 1 ) {
// 	  //std::cout << "sources->Ecut() " << sources->Ecut()<< std::endl;
// 	  _fEnergy = pow(sources->Emin()/double(Amax),toto)*pow(sources->Ecut(),1-toto) ;
// 	} else {
// 	  _fEnergy = pow ( pow(sources->Emin()/double(Amax),1-alpha) + toto * 
// 			   (pow(sources->Ecut(),1-alpha)-pow(sources->Emin()/double(Amax),1-alpha)) , 1/(1-alpha) ) ;
// 	  cout<<_fEnergy*invEeV<<" new "<< invEeV* RandPowerLaw(-alpha, sources->Emin()*double(Amax),sources->Ecut()*Amax,toto)<<endl;
// 	}
	_fEnergy=RandPowerLaw(-alpha, sources->Emin(),lEmax,toto);
      }
    } else throw TCrpErr("Unknown source type in TNucleus constructor" );
//    if(sources->RigidityFlag()) _fEnergy *= _fChargeNumber; //In case of rigidity multiply by charge
  
  } while( sources->RigidityFlag()&& (_fEnergy < sources->Emin() || _fEnergy/double(_fChargeNumber) > sources->Ecut()));


  // 2 ) Initial momentum
  if ( aUniv->Type() == UNIVERSE_ENV1D) {

    _fMomentum.setX( - _fEnergy* inv_c_light ) ;
    //std::cout<<"initial momentum 1d"<<std::endl;
  } else { // 3D random direction
    _fMomentum.setRThetaPhi( _fEnergy* inv_c_light, acos(1-2*RandFlat::shoot()), twopi * RandFlat::shoot() ) ;

    /*Testing
    std::cout<< "momentum set along X axis" << std::endl;
    _fMomentum.setX( _fEnergy* inv_c_light ) ;
    _fMomentum.setY(0 ) ;
    _fMomentum.setZ(0 ) ;
    Testing*/ 

  }
  _fInitPosition = _fPosition ;
  _fInitMomentum = _fMomentum ;
  interactt+=3;

    
  // 3 ) Initial timesteps and distance to observers
  _fDistObs = 1 * Gpc ;
  if (aUniv->Observers()->Type() == OBSERVER_SMALLSPHERE ) {
    TObservers* lObs = aUniv->Observers() ;
    _fWasInsideSpheres.assign(aUniv->Observers()->Nb(),false) ;
    for (int j=0; j<lObs->Nb(); j++) {
      double lDistObs = (_fPosition - lObs->Positions(j)).mag() - lObs->Radius() ;
      if (lDistObs > 0) _fDistObs = min(_fDistObs,lDistObs) ;
    }
  }
  if (aUniv->Observers()->Type() == OBSERVER_LARGESPHERE )
    _fIsInsideSpheres.assign(aUniv->Observers()->Nb(),1) ;


   _fTimeStep = aUniv->IntegratorMinTimeStep() ;
   _fNextTimeStep = _fTimeStep; 
 
  // Initial redshift
  if (aUniv->Type() != UNIVERSE_ENV1D || !_fpUniverse->InteractionData()->RedshiftFlag()) {
    _fRedshift = 0 ;
  } else{
    this->ComputeRedshift() ;
    _fInitRedshift=_fRedshift;
  }
  
  TInteractionData *lpInt = _fpUniverse->InteractionData() ;

  //  std::cout<<" lpInt->RedshiftFlag() " <<  lpInt->RedshiftFlag() << " _fRedshift " << _fRedshift << std::endl;  
  // #ifdef HAVE_TROOT_H


  ////////////////////////////////////////////////////////////////////
  // Routines to plot photon densities (begin)
    // please, comment in when needed
  /*
  CMB CMB_B;
  TGraph * CMBGraph = RootPlotSpectrum(CMB_B,
 				       4.e-19,
 				       12.4e-9,
				       500,
				       0.,
				       0.,
				       0.,
				       0.);
  TCanvas* CMBCan = new TCanvas("CMBCan","CMBCan",1);
  CMBGraph->Draw("A*");
  CMBCan->SaveAs("CMBCan.root");
  */
  
  //IR IR_B;
  /*TVariableInfrared* IR_B = new TVariableInfrared("/home/nils/software/CRPropa/trunk/examples/3dSimTest/plotMFP.xml");
  TGraph * IRGraph = RootPlotSpectrum((*IR_B),
 				       4.e-19,
 				       12.4e-9,
				       500,
				       10.,
				       -1.,
				       -1.,
				       0.);
  TCanvas* IRCan = new TCanvas("IRCan","IRCan",1);
  IRGraph->Draw("A*");
  IRCan->SaveAs("IRCan.root");
  */
  // Routines to plot photon densities (end)
  ////////////////////////////////////////////////////////////////////



  ////////////////////////////////////////////////////////////////////////////////////
  // Create & store root plots for single channels from photo disintegration tables. (begin)
  // please, comment in when needed
  //TabulatedTALYSCrossSection::GetInstance()->RootPlotTotalY(26056,26055);
  //TabulatedTALYSCrossSection::GetInstance()->ASCIIOutExclY(10025,212101);
  //TabulatedTALYSAveragedCrossSection::GetInstance()->RootPlotTotalY(26056,26055);
  //TabulatedTALYSAveragedCrossSection::GetInstance()->RootPlotTotalY(26056,26055);
  //TabulatedTALYSAveragedCrossSection::GetInstance()->RootPlotExclY(2006,100000);
  //TabulatedTALYSMeanFreePath::GetInstance()->RootPlotTotalY(26056,26055);
  //TabulatedTALYSMeanFreePath::GetInstance()->RootPlotExclY(26056,110000);
  //exit(-1);
  // Create & store root plots for single channels from photo disintegration tables. (end)
  ////////////////////////////////////////////////////////////////////////////////////



  ////////////////////////////////////////////////////////////////////////////////////
  // Create ASCII files for single channels from photo disintegration tables. (begin)
  // please, comment in when needed
  //TabulatedTALYSCrossSection::GetInstance()->ASCIIOutPlotTotalY(26056,26055);  
  //TabulatedTALYSCrossSection::GetInstance()->ASCIIOutExclY(26056,110000);
  //TabulatedTALYSAveragedCrossSection::GetInstance()->ASCIIOutPlotTotalY(19041,15029);  
  //TabulatedTALYSAveragedCrossSection::GetInstance()->ASCIIOutExclY(26056,100000);
  //TabulatedTALYSMeanFreePath::GetInstance()->ASCIIOutPlotTotalY(26056,26055);  
  //TabulatedTALYSMeanFreePath::GetInstance()->ASCIIOutExclY(26056,100000);
  //exit(-1);
  // Create ASCII files for single channels from photo disintegration tables. (end)
  ////////////////////////////////////////////////////////////////////////////////////
  


#ifdef PLOT_MEAN_FREE_PATH_ROOT 
  /////////////////////////////////////////////////////////////////////////////////////
  // P L O T     R O U T I N E S     F O R     P D     M E A N     F R E E     P A T H / XS (begin)
  // loop over all channels included 
  //CASE A: This plots the mean free path from the mfp tables 

  std::cout<<"\n\n***Compiled with -DPLOT_MEAN_FREE_PATH_ROOT to create plots for photo dosintegration.:"<<std::endl;
  std::cout<<"\n***Select the output:"<<std::endl;
  std::cout<<"\t1. mean free path tables (fast),"<<std::endl;
  std::cout<<"\t2. mean free path from averaged cross section (slow),"<<std::endl;
  std::cout<<"\t3. cross section tables."<<std::endl;
  std::cout<<"The default is 1."<<std::endl;
  int inputCase;
  std::cin>>inputCase;
  if(inputCase<1 || inputCase>3) inputCase=1;
  
  double gammaMin=0.;
  double gammaMax=0.;
  int NPoints=0.;
  double redshift=0.;
  TVector3D Position; Position.setX(0.);Position.setY(0.);Position.setZ(0.);
  double epsilon0=0., epsilon_min=0., epsilon_max=0.;

  if(inputCase==1 || inputCase==2) {
    std::cout<<"Give minimum log10 gamma factor (e.g. 6):"<<std::endl;
    std::cin>>gammaMin;
    std::cout<<"Give maximum log10 gamma factor (e.g. 14):"<<std::endl;
    std::cin>>gammaMax;
    std::cout<<"Give number of points (e.g. 200):"<<std::endl;
    std::cin>>NPoints;
    std::cout<<"Give redshift:"<<std::endl;
    std::cin>>redshift;
  }


  if(inputCase==1){
    TALYSMeanFreePathTabulated TALYSMeanFreePathTabulated;
    std::vector<PhotonBackground*> PhotonBackgrounds_Vec; 
    std::vector<TabulatedTALYSY*> TabData_Vec; 
    int input=0;
    std::cout<<"Choose the photon field:"<<std::endl;
    std::cout<<"1. CMB"<<std::endl;
    std::cout<<"2. IRB"<<std::endl;
    std::cout<<"3. CMB and IRB"<<std::endl;
    std::cin>>input;
    if(input==1 || input ==3) TabData_Vec.push_back(TabulatedTALYSMeanFreePath_CMB::GetInstance());
    if(input==2 || input ==3) TabData_Vec.push_back(TabulatedTALYSMeanFreePath_IRB::GetInstance());
    TALYSMeanFreePathTabulated.RootFileWithAllMFPs(PhotonBackgrounds_Vec, 
						   TabData_Vec, 
						   Position, 
						   gammaMin, 
						   gammaMax, 
						   NPoints, 
						   0.);
    exit(-1);
  }
  
  if(inputCase==2){


    std::cout<<"give epsilon_min (e.g. 4.e-19):"<<std::endl;
    std::cin>>epsilon_min;
    std::cout<<"give epsilon_max (e.g. 12.4e-9):"<<std::endl;
    std::cin>>epsilon_max;
    std::cout<<"give epsilon0 (e.g. kT=8.617e-14*2.725=2.348e-13):"<<std::endl;
    std::cin>>epsilon0;
    
    TALYSMeanFreePathAvrgd TALYSMeanFreePathAvrgd(epsilon0, epsilon_min, epsilon_max);
    std::vector<PhotonBackground*> PhotonBackgrounds_Vec;

    int input=0;
    std::cout<<"Choose the photon field:"<<std::endl;
    std::cout<<"1. CMB"<<std::endl;
    std::cout<<"2. IRB"<<std::endl;
    std::cout<<"3. CMB and IRB"<<std::endl;
    std::cin>>input;
    
    if(input==1 || input ==3) PhotonBackgrounds_Vec.push_back(new CMB());
    
    if(input==2 || input ==3){      
      std::cout<<"\nWhat IRB you wanna use? \n1. const IRO(Primack)\n2. variable infrared\n>>"<<std::endl;
      int VarOrConsIRB=1;
      std::cin>>VarOrConsIRB;
      if(VarOrConsIRB!=2){
	PhotonBackgrounds_Vec.push_back(new IR());
      }
      else{
	std::cout<<"Pleas give the path to the xml configuration file where the variable IRO is defined:"<<std::endl;
	std::string FileName;
	std::cin>>FileName;
	
	PhotonBackgrounds_Vec.push_back(new TVariableInfrared(FileName.c_str(), true));  
	
	std::cout<<"Position x/Mpc:"<<std::endl;
	double bufferX=0.;
	std::cin>>bufferX;
	Position.setX(bufferX*Mpc);
	std::cout<<"Position y/Mpc:"<<std::endl;
	double bufferY=0.;
	std::cin>>bufferY;
	Position.setY(bufferY*Mpc);
	std::cout<<"Position z/Mpc:"<<std::endl;
	double bufferZ=0.;
	std::cin>>bufferZ;
	Position.setZ(bufferZ*Mpc);
	
	std::cout<<"Position is:"
		 <<" x="<<Position.getX()/Mpc
		 <<" y="<<Position.getY()/Mpc
		 <<" z="<<Position.getZ()/Mpc
		 <<std::endl; 
      }
    }
    
    std::vector<TabulatedTALYSY*> TabData_Vec; 
    TabData_Vec.push_back(TabulatedTALYSAveragedCrossSection::GetInstance());
    TALYSMeanFreePathAvrgd.SetNIntegratorCalls(15);
    TALYSMeanFreePathAvrgd.RootFileWithAllMFPs(PhotonBackgrounds_Vec, 
					       TabData_Vec,
					       Position, 
					       gammaMin, 
					       gammaMax, 
					       NPoints, 
					       0.);
    exit(-1);
  }
  
  if(inputCase==3){
    //CASE C: This plots the cross sections 
    double Emin=TabulatedTALYSCrossSection::GetInstance()->GetMinimumXValue()-5.e-4;
    double Emax=TabulatedTALYSCrossSection::GetInstance()->GetMaximumXValue()+10;
    TabulatedTALYSCrossSection::GetInstance()->RootFileWithAllReactions(Emin, Emax, 500, 0);
    exit(-1);
  }
  // P L O T     R O U T I N E S     F O R     P D     M E A N     F R E E     P A T H (end)
  /////////////////////////////////////////////////////////////////////////////////////
#endif

} // End of constructor




/*************************************************************/
QUEUE<TParticle*>* TNucleus::Propagate(TBasicParam* aBasic) {

    TInteractionData *lpInt = _fpUniverse->InteractionData() ;
    double TestDistanceSum=0.;

    /* This code prints the mfp for all interactions into the file mfpPPOut.txt.  */
    //#define mfpOut
#ifdef mfpOut
    fstream mfpOutFile ("mfpPPOut.txt", fstream::out);
   if( mfpOutFile.fail()) throw TCrpErr("mfpOutfile Fail!");
   _fRedshift=0;
#endif
 

    // CreatePionTable();
  if (aBasic->RecordMode() == "Full Trajectories" && _fpUniverse->Type() == UNIVERSE_LARGESCALE)
    aBasic->OutputData()->Add3DTraj(-1,-1 * Mpc/ c_light ,_fInitPosition,_fInitMomentum,-1 *EeV, -1) ;
  //Allows to separate individual trajectories as well as to give the true initial position

  do {

    if (aBasic->RecordMode() == "Full Trajectories") this->Write(aBasic) ;
    if (aBasic->RecordMode() == "Events") this->CheckDetection(aBasic) ;

    if ((_fpUniverse->Type() == UNIVERSE_ENV1D && lpInt->RedshiftFlag())) {
      this->ComputeRedshift() ;
      if(lpInt->RedshiftEnergyLossFlag()) this->RedshiftLoss() ;
    }
    
    
    if (lpInt->Type() == INTERACTION_NO) {
      //Nothing to be done for no interaction. 
      if ( _fpUniverse->Type() == UNIVERSE_ENV1D ) _fDistObs=_fPosition.x() ;
      _fNextTimeStep=min(_fpUniverse->InteractionData()->InteractionTimeStep(),.9 * _fDistObs* inv_c_light);
      this->Deflec() ;   // Propagates the particle through *aUniv
    } else if (
	       
	       lpInt->Type() == INTERACTION_BASIC 
	       &&(_fMassNumber == 1)//TODO : add this to the consistency
	       ) {
      throw TCrpErr("Interactions F77-proton no longer supported. Use sophia interaction type instead"); 
      if(_fNextTimeStep > 0.)this->Deflec() ;   // Propagates the particle through *aUniv
      if (lpInt->PairProdFlag()) this->PairProduction( _fNextTimeStep) ;
      std::cout << " INTERACTION_BASIC UEBERARBEITEN!!!" << std::endl; 
      if (lpInt->PionProdFlag()) this->PionInteraction() ;
      //      this->NucleusDecay();
      // Here the InteractionStepSize is set like in CRPropa 1.3
 
   } else
    
  if (lpInt->Type() == INTERACTION_SOPHIA) {
      
      //Get the mean free path of all reactions. Order:
      //1. Nucleus Decay
      //2. Photodisintegration
      //3. PionProduction on proton (CMB)
      //4. PionProduction on neutron (CMB)
      //5. PionProduction on proton(IR)
      //6. PionProduction on neutron(IR)
      //7. Proton-Proton on target gas
      //TODO check all units
     
    //These are the energies before Pairproduction and redshift energz losses
    double EnergyAtStartPoint=_fEnergy;
    double RedshiftAtStartPoint=_fRedshift;
    
    std::vector<double> MFPOfAllReactions;
    double  PDAllChannelMFP;
    
    
    /*TestCode */  
#ifdef mfpOut
      double OverallMFP=0.;
      mfpOutFile << "#Energy Gamma Decay  PD PionProtonCMB PionNeutronCMB PionProtonIRO PionNeutronIRO ProtonProtonOnTargetGas MfP PairPRodLossLength @z="<< _fRedshift << std::endl; 
      
      for( _fEnergy=.01 *EeV; _fEnergy < 1e5 * EeV; _fEnergy *= pow(10., 1/25.) ){
	OverallMFP=0.;
	MFPOfAllReactions.clear();
          MFPOfAllReactions.reserve(7);
	std::cout << " _fEnergy : " << _fEnergy << std::endl;
#endif
	
	TVector3D Position; Position.setX(0.);Position.setY(0.);Position.setZ(0.);

	////////////////////////////////////////////////////////////////////////
	//1. Decay (begin)

#ifdef DEBUG_OUTPUT
	static int lFlag=1;
	if(lFlag){
	  std::cout << "DecaYFlag " << lpInt->DecayFlag() << std::endl;
	  lFlag=0;
	}
	std::cout << "DecayGammaTau()" ;
#endif
	if(lpInt->DecayFlag()) {
	  MFPOfAllReactions.push_back(DecayGammaTau());
	} else {
	  //dummy entry to guarantee a vector with 6 components
	  MFPOfAllReactions.push_back(std::numeric_limits<double>::max());	
	}

#ifdef DEBUG_OUTPUT
	std::cout << "done" << std::endl;
#endif
	//1. Decay (end)
	////////////////////////////////////////////////////////////////////////


	////////////////////////////////////////////////////////////////////////
	//2.Photodisintegration (begin)
#ifdef DEBUG_OUTPUT
	std::cout << "PD" ;
#endif
	if(_fMassNumber >= 2 && lpInt->PhotodisintegrationFlag()){
	  PDAllChannelMFP = 1./lpInt->GetPDTimeStep(_fChargeNumber*1000 + _fMassNumber,
						    _fPosition,
						    _fEnergy*(1.+_fRedshift)/1000.,  //GeV
						    _fMass*c_squared/1000.,             //GeV
						    _fRedshift );
	  MFPOfAllReactions.push_back(PDAllChannelMFP);
	  //std::cout<<"PDAllChannelMFP in Prop() :"<<PDAllChannelMFP<<std::endl;
#ifdef DEBUG_OUTPUT
	std::cout << "done" << std::endl;
#endif
	}else{
	  //dummy entry to guarantee a vector with 6 components
	  MFPOfAllReactions.push_back(std::numeric_limits<double>::max());	
	}
      	//2.Photodisintegration (end)
	////////////////////////////////////////////////////////////////////////
	


	////////////////////////////////////////////////////////////////////////
	//3./4. Pion production on CMB (p/n) (begin)
#ifdef DEBUG_OUTPUT
	std::cout << "PionProton /Neutron CMB " << std::endl;
#endif
	if (lpInt->PionProdFlag()){
	  MFPOfAllReactions.push_back(PionProtonMFP( 1));
	  MFPOfAllReactions.push_back(PionNeutronMFP( 1));
	}else{
	  //dummy entry to guarantee a vector with 6 components
	  MFPOfAllReactions.push_back(std::numeric_limits<double>::max());
	  MFPOfAllReactions.push_back(std::numeric_limits<double>::max());
	}
#ifdef DEBUG_OUTPUT
	std::cout << "done" << std::endl;
#endif
	//3./4. Pion production on CMB (p/n) (end)
	////////////////////////////////////////////////////////////////////////
	
	
	////////////////////////////////////////////////////////////////////////
	//5./6. Pion production on IRO (p/n) (begin)
#ifdef DEBUG_OUTPUT
	std::cout << "Pion Proton/Neutron done" << std::endl;
#endif
	if (lpInt->IRPionProdFlag() && lpInt->PionProdFlag()){
	  if (
	      (_fpUniverse->Infrared()->Type() == SHELL || _fpUniverse->Infrared()->Type() == THREED) 
	      && 
	      _fMassNumber == 1
	      ) {
	    MFPOfAllReactions.push_back(PionProtonMFP( 3));
	    MFPOfAllReactions.push_back(PionNeutronMFP( 3));
	  }
	  else {
	    MFPOfAllReactions.push_back(PionProtonMFP( 2));  
	    MFPOfAllReactions.push_back(PionNeutronMFP( 2));  
	  }
	} else{
	  //dummy entry to guarantee a vector with 6 components
	  MFPOfAllReactions.push_back(std::numeric_limits<double>::max());
	  MFPOfAllReactions.push_back(std::numeric_limits<double>::max());
	}
#ifdef DEBUG_OUTPUT
	std::cout << "done" << std::endl;
#endif
	//5./6. Pion production on IRO (p/n) (end)
	////////////////////////////////////////////////////////////////////////

	
	////////////////////////////////////////////////////////////////////////
	//7. Proton-Proton on target gas (begin)
	if (lpInt->PProdFlag() && _fMassNumber == 1){ 
	  MFPOfAllReactions.push_back(SibyllInteraction());
	}
	else {
	  //dummy entry to guarantee a vector with 6 components
	  MFPOfAllReactions.push_back(std::numeric_limits<double>::max());
	}
	//7. Proton-Proton on target gas (end)
	////////////////////////////////////////////////////////////////////////
	

      //Inversly sum up to one overall mean free path.
      double OverallMFP=0.;
      /*TESTCODE  */
#ifdef mfpOut
      try{
	mfpOutFile << _fEnergy /EeV << "  " <<  _fEnergy / this->Mass() / c_squared;
#endif
	
	for(int i=0; i<MFPOfAllReactions.size(); i++){
	  /*TESTCODE   */
#ifdef mfpOut
	  mfpOutFile << "  " << MFPOfAllReactions[i]/(Mpc * inv_c_light);
#endif
	  if(MFPOfAllReactions[i]!=0.){
	    OverallMFP+=1./MFPOfAllReactions[i];
	    
	  }else{
	    OverallMFP+=1./1.e99;
	    std::cout<<"Error: MFPOfAllReactions["<<i<<"]=0"<<std::endl;  
	    throw TCrpErr("ERROR: MFPOfAllReactions[i]=0");  
	    }
	  
	}
	OverallMFP=1/OverallMFP;
#ifdef mfpOut
	mfpOutFile << "  " <<  OverallMFP  /(Mpc * inv_c_light) << "  " << PairProdLossLength()/(Mpc * inv_c_light) <<std::endl; 
      } catch (exception& e){
	std::cout << " e.what() " <<  e.what() << std::endl;
	std::cout << " _fEnergy " << _fEnergy << std::endl;
      }
      }
      
      mfpOutFile.close();
      throw TCrpErr("STOP");
#endif

      
      //Find the distance to the next reaction via MC (use inverse decay law.)
      //F.Y.I.: In what follows now we use the max time step from SophiaInteraction (->InteractionData()->InteractionTimeStep()). And min time step for Runge-Kutta Integrator (aUniv->IntegratorMinTimeStep()).
      double p                  = RandFlat::shoot();
      double DistanceToReaction = -1.*OverallMFP*log(p);
      bool PerformReaction=false;
      //#define PairProd_EPS 0.01
 
     static const double PairProd_EPS= _fpUniverse->InteractionData()->GetPPEps();
     if ( _fpUniverse->Type() == UNIVERSE_ENV1D ) _fDistObs=_fPosition.x() ;
     double hmaxstor=min(_fpUniverse->InteractionData()->InteractionTimeStep(),.9 * _fDistObs* inv_c_light);
     if (_fpUniverse->InteractionData()->PairProdFlag())  hmaxstor=min( PairProd_EPS * PairProdLossLength(), hmaxstor);
     _fNextTimeStep=min(DistanceToReaction, hmaxstor);
     if( hmaxstor >= DistanceToReaction) PerformReaction=true; 
     
     //    if(_fNextTimeStep < numeric_limits<double>::min() ) {
     if(_fNextTimeStep ==0 ){
       
       std::cout << "PP Losslength : " <<  PairProd_EPS * PairProdLossLength() /(Mpc * inv_c_light)  << std::endl;
       std::cout << ".9 * _fDistObs* inv_c_light  " << .9 * _fDistObs* inv_c_light   /(Mpc * inv_c_light) << std::endl; 
       std::cout << " DistanceToReaction " << DistanceToReaction  /(Mpc * inv_c_light) << std::endl; 
       std::cout << " OverallMFP " << OverallMFP   /(Mpc * inv_c_light) << std::endl; 
       
       std::cerr<< " WARNING : Nucleus : _fNextTimeStep < numeric_limits<double>::min() " << std::endl; 
       std::cerr<< " _fNextTimeStep =" << _fNextTimeStep/(Mpc * inv_c_light) << "  numeric_limits<double>::min()=" <<numeric_limits<double>::min()<<std::endl;
       std::cout<< "PairProd_EPS " << PairProd_EPS << std::endl;
       std::cout<<"p="<<p<<" DistanceToReaction="<<DistanceToReaction/(Mpc * inv_c_light)<<std::endl;
       std::cout<<"OverallMFP="<< OverallMFP<<std::endl;
       for(std::vector<double>::iterator i =MFPOfAllReactions.begin(); i < MFPOfAllReactions.end(); i++) std::cerr<< *i /(Mpc * inv_c_light) << std::endl;
       std::cout << "_fRedshift : " << _fRedshift << std::endl;
     }
     
     
#ifdef DEBUG_OUTPUT
     std::cout << "Deflec";
#endif
     this->Deflec() ;
     if (_fpUniverse->InteractionData()->PairProdFlag()) this->PairProduction( _fTimeStep) ;
     // Propagates the particle through *aUniv
     
     TestDistanceSum+=_fNextTimeStep;
     
     if(_fMassNumber==56 && _fChargeNumber==26){ TestDistanceSum_Excl+=_fNextTimeStep;}

#ifdef DEBUG_OUTPUT
     std::cout << " done" << std::endl;
#endif
     //Choose a channel
     int k;
     //    PerformReaction=false;
     
     

     if(PerformReaction==true){
       

       //if(DistanceToReaction < _fDistObs* inv_c_light
       //	 &&
       //	 DistanceToReaction < _fpUniverse->InteractionData()->InteractionTimeStep()
       //	 ){
      
      double probSum   = 0.;
      double RanNumber = RandFlat::shoot();
      
      for(k=0; k<MFPOfAllReactions.size(); k++){
	probSum += (OverallMFP/MFPOfAllReactions[k]); 
	//	  std::cout<<"k="<<k<<"\tprobSum="<<probSum<<std::endl;
	if(probSum>=RanNumber) break;
      }
	
       
#ifdef DEBUG_OUTPUT
      std::cout << "Chosen reaction: " << k << std::endl;
#endif
      bool NilsPlotFlag=false;
      if(_fMassNumber==56 && _fChargeNumber==26 && k==1){
	NilsPlotFlag=true;
	//if(_fRedshift<=0.5) NilsPlotFlag=true;
      }
      
      double storE=_fEnergy; 
      //Write the reactions products to particle queue.
      //Decay
      if(k==0){
	GetNuclearDecayProducts();
	//NeutronDecay();
      }
      //Photodisintegration;    
      else if(k==1){	  
          
	//MassChargePathTH2F->Fill(_fChargeNumber+.5, _fMassNumber+.5);	  
	this->TALYSNucleiPhotoDisintegration(aBasic);
	//MassChargePathTH2F->Fill(_fChargeNumber+.5, _fMassNumber+.5);	  

#ifdef PDSEC
	// New electromagnetic cascade.                                  
	//Just a very quick hack to have PD secondary gammas
	double lCharge = 0;
	TVector3D lGammaMomentum = _fMomentum ;
	lGammaMomentum.setMag(10 * MeV / c_light * _fEnergy/ this->Mass()/ c_squared) ;
	if (lpInt->SecPhotonFlag())
	  if (_fpUniverse->Type() == UNIVERSE_ENV1D) {
	    PARTICLE	lPart = PHOTON ;
	    _fpList1DPhotons->AddPhotons(_fPosition.x(), lGammaMomentum.mag()*c_light, lPart) ;
	  } else { // Universe Largescale
	    string lPart ;
	    if (lCharge == 0) {
	      lPart = "pd_gamma";
	    } else throw TCrpErr("Charge error!") ;
	    fParticleQueue.push_back(new TPhoton(_fPosition, lGammaMomentum, lCharge, _fInitPosition,
						 _fInitMomentum.getR()*c_light, _fpUniverse, _fInitType, lPart, _fTime));
	  }
#endif
	}
	//Photopionproduction of a proton on the CMB
	else if(k==2){
	  this->SophiaPionProd( 1, 1, _fEnergy / (double) _fMassNumber) ;
	}
	//Photopionproduction of a neutron on the CMB
	else if(k==3){
	  this->SophiaPionProd( 1, 0, _fEnergy / (double) _fMassNumber) ;
	}
	//Photopionproduction of a proton on the IRO

	else if(k==4){
	  if ((_fpUniverse->Infrared()->Type() != SHELL && _fpUniverse->Infrared()->Type() != THREED) || _fMassNumber > 1){
	    this->SophiaPionProd( 2, 1, _fEnergy / (double) _fMassNumber) ;
	  }
	  else{ 
	    this->SophiaPionProd( 3, 1, _fEnergy / (double) _fMassNumber) ;	    
	  }
	}
	//Photopionproduction of a neutron on the IRO
	else if(k==5){
	  if ((_fpUniverse->Infrared()->Type() != SHELL &&  _fpUniverse->Infrared()->Type() != THREED) || _fMassNumber > 1)
	    this->SophiaPionProd( 2, 0, _fEnergy / (double) _fMassNumber) ;
	  else 
	    this->SophiaPionProd( 3, 0, _fEnergy / (double) _fMassNumber) ;	    
	}  
	// Proton-proton on gas
	else if (k==6) {
	  this->SibyllPProd();
	}
#ifdef DEBUG
      if(_fEnergy > storE +1.e3 ){//Test if the particle gains energy
	    std::cout<<"k="<<k<<std::endl;
	    std::cout<<"_fEnergy="<<_fEnergy<<std::endl;
	    std::cout<<"_fMassNumber="<<_fMassNumber<<std::endl;
	  throw  TCrpErr("PairProd: _fEnergy > storE + 1 " ); 
	}
         
         
      //Test for inconsistency between Charge Number and Charge
	if( (_fChargeNumber - _fCharge/eplus) > 1e-3)std::cerr << "  eplus *_fChargeNumber - _fCharge < .1 : Channel " << k << " _fChargeNumber " << _fChargeNumber << "_fCharge / eplus " << _fCharge / eplus << std::endl
#endif
	if(NilsPlotFlag==true){
#ifdef NILS
	  //MassChargeVSChannelTH2F->Fill(log10(_fEnergy/(_fMass*c_squared)), k);
	  //MFPvsEnTH2F->Fill(log10(_fEnergy/(_fMass*c_squared)), log10(DistanceToReaction/(Mpc * inv_c_light)));
	  //std::cout<<"Test:"<<log10(_fEnergy/(_fMass*c_squared))<<" "<<log10(TestDistanceSum/(Mpc * inv_c_light))<<std::endl;
	  if(_fChargeNumber==25 && _fMassNumber==55) {
	    MFPvsEn_ExclTH2F->Fill(log10(_fEnergy/(_fMass*c_squared)), log10(TestDistanceSum_Excl/(Mpc * inv_c_light)));
	    TestDistanceSum_Excl=0.; 
	    //std::cout<<"case 1n"<<std::endl;
	  }//else{std::cout<<"case !=1n"<<std::endl;}
	  MFPvsEnTH2F->Fill(log10(_fEnergy/(_fMass*c_squared)), log10(TestDistanceSum/(Mpc * inv_c_light)));
	  //std::cout<<"Fill now"<<std::endl;
#endif
	}
	TestDistanceSum=0.;
    }
  } else throw TCrpErr("Proton interactions not implemented" );  
#ifdef NILS
    //MassChargePathTH2F->Fill(_fChargeNumber+.5, _fMassNumber+.5); 	  
#endif

 #ifdef DEBUG_OUTPUT
    std::cout << "!!!End of propagation step!!! " << std::endl;
#endif
  } while (this->CheckEndPropa(aBasic) == false);

  //#define SPEED_TEST
#ifdef SPEED_TEST
  static fstream speedtst ("speed.dat", fstream::out);
  static clock_t clockTime=0;
  static clock_t oldTime=0;
  static struct timeval tv0, tv1;

  if( speedtst.fail()) throw TCrpErr("speed.dat Fail!");
  
  oldTime=clockTime;
  clockTime=clock();
  tv0.tv_sec=tv1.tv_sec;
  tv0.tv_usec=tv1.tv_usec;
  gettimeofday(&tv1, NULL);

  speedtst << _fInitMomentum.getR()*c_light / EeV << "  " << (tv1.tv_sec-tv0.tv_sec)*1e6+(tv1.tv_usec-tv0.tv_usec) << "  " << (double)(tv1.tv_sec*1e6 + tv1.tv_usec) << "  " << (clockTime-oldTime)/(double)CLOCKS_PER_SEC <<  std::endl;
  
#endif 

  return &fParticleQueue;
}



double TNucleus::PairProdLossLength(){
  TInteractionData *lpInt = _fpUniverse->InteractionData() ;


  if(abs(_fCharge)){
    double norm_p = _fEnergy * invEeV * (1.+_fRedshift); // f77 units -- REDSHIFT ADDED!
  norm_p /= this->Mass() * c_squared / proton_mass_c2;
#define USE_PAIRPRODRATE_EXTENDED
#ifdef USE_PAIRPRODRATE_EXTENDED
    int index_E = min((int)floor(log10(norm_p/lpInt->E_part(0)) /lpInt->dEtabbin()) , 69);
#else
    int index_E = min((int)floor(log10(norm_p/lpInt->E_part(0)) /lpInt->dEtabbin()) , 79) ;
#endif

      if (index_E >= 0 ) { // check out indexes.
	double lossLength;
	if(index_E < 69){
	  lossLength=
	pow(1.+_fRedshift,3.)* proton_mass_c2 / _fMass / c_squared * _fChargeNumber * _fChargeNumber *
	(
	 (
	  (norm_p-lpInt->E_part(index_E)
	   )
	  *
	  (lpInt->L_pair(index_E+1)-lpInt->L_pair(index_E))
	  /
	  (lpInt->E_part(index_E+1)-lpInt->E_part(index_E))
	  )
	 +lpInt->L_pair(index_E)
	 ) ; // REDSHIFT ADDED
	}else {
	  lossLength=
	    pow(1.+_fRedshift,3.)* proton_mass_c2 / _fMass / c_squared * _fChargeNumber * _fChargeNumber * lpInt->L_pair(index_E) * pow(norm_p / lpInt->E_part(index_E), -.6);// Linear interpolation: log(Losslength)=14430 Mpc * (Energy/ 1e22 eV)^.6
	}

      if(index_E == 70){
	std::cerr<< "PairProdLossLength: WARNING: index_E max!" << std::endl;
	std::cerr << "index_E " << index_E << std::endl;
	std::cerr<< "_fRedshift=" << _fRedshift << "\t_fEnergy=" << _fEnergy << std::endl; 
	std::cerr << " lpInt->E_part(0)" <<  lpInt->E_part(0) << "\tlpInt->dEtabbin() " << lpInt->dEtabbin() << std::endl;
	std::cerr << "this->Mass() * c_squared / proton_mass_c2 " << this->Mass() * c_squared / proton_mass_c2 << std::endl;
	exit(-1);
      }
      if( lossLength < 0){
	std::cout << "PP LL < 0 " << lossLength << std::endl;
	std::cout << "_fCharge " << _fChargeNumber << std::endl;
	std::cout << "_fMassNumber " << _fMassNumber << std::endl;
	std::cout << "1rst Line " << pow(1.+_fRedshift,3.)* proton_mass_c2 / _fMass / c_squared * _fChargeNumber * _fChargeNumber << std::endl;
	std::cout << "2nd Line " << (norm_p-lpInt->E_part(index_E)) << std::endl;
	std::cout << "3rd " << (lpInt->L_pair(index_E+1)-lpInt->L_pair(index_E)) << std::endl;
	std::cout << "4th " << (lpInt->E_part(index_E+1)-lpInt->E_part(index_E)) << std::endl;
	std::cout << "5th " << lpInt->L_pair(index_E) << std::endl;
	std::cout << "index_E " << index_E << std::endl;
	std::cout << "norm_p " << norm_p << std::endl; 
	std::cout << "_fEnergy * invEeV " << _fEnergy * invEeV << std::endl;
	std::cout << "lpInt->E_part(index_E) " << lpInt->E_part(index_E) << std::endl;
	exit(-1);
      }
    
      return 1/lossLength * Mpc / c_light;
      
    }


  }
  

  return 1e99;
}

void TNucleus::PairProduction(double lTimeStep) {
// pair production from f77 code  
  if( lTimeStep == 0) return; 
 TInteractionData *lpInt = _fpUniverse->InteractionData() ;
 static int testflag=0; 
 if(testflag==0){
  testflag=1;
  std::cout<< "\t->PairProduction called" << std::endl; 
  //std::cout<< "Achtung modified!!! " << std::endl; 
 }

 
  if (abs(_fCharge)) {
    double dfin = (lTimeStep * c_light) * invMpc ; // f77 units
    double norm_p = _fEnergy * invEeV * (1.+_fRedshift); // f77 units -- REDSHIFT ADDED!
    norm_p /= this->Mass() * c_squared / proton_mass_c2;
#define USE_PAIRPRODRATE_EXTENDED
#ifdef USE_PAIRPRODRATE_EXTENDED
    int index_E = min((int)floor(log10(norm_p/lpInt->E_part(0)) /lpInt->dEtabbin()) , 69);
#else
    int index_E = min((int)floor(log10(norm_p/lpInt->E_part(0)) /lpInt->dEtabbin()) , 79) ;
#endif

    if (index_E >= 0 ) { // check out indexes...
      double loss_pair;
      if( index_E < 69){
      loss_pair =
	pow(1.+_fRedshift,3.)*dfin* proton_mass_c2 / _fMass / c_squared * _fChargeNumber * _fChargeNumber *
	(
	 (
	  (norm_p-lpInt->E_part(index_E)
	   )
	  *
	  (lpInt->L_pair(index_E+1)-lpInt->L_pair(index_E))
	  /
	  (lpInt->E_part(index_E+1)-lpInt->E_part(index_E))
	  )
	 +lpInt->L_pair(index_E)
	 ) ; // REDSHIFT ADDED
      } else {
	loss_pair=	pow(1.+_fRedshift,3.)*dfin* proton_mass_c2 / _fMass / c_squared * _fChargeNumber * _fChargeNumber *lpInt->L_pair(index_E) * pow(norm_p / lpInt->E_part(index_E), -.6);// Linear interpolation: log(Losslength)=14430 Mpc * (Energy/ 1e22 eV)^.6 
      }
      // for cmb-like redshift evolution of the background, the energy loss coefficient 
      // is dE/dt = b(E,z) = (1+z)^2 b(E(1+z)) [Bugaev]
      // This implies (1/E dE/dt)(E,z) = (1+z)^3 (1/E dE/dt)(E(1+z))

      if (loss_pair >= 1) throw TCrpErr("Too large stepsize for pair production") ;
      _fMomentum *= 1-loss_pair ;
      if(loss_pair == 0){
	std::cerr << " loss_pair == " <<loss_pair << "   _fEnergy=" << _fEnergy << std::endl;
	std::cerr << " pow(1.+_fRedshift,3.)*dfin* proton_mass_c2 / _fMass / c_squared * _fChargeNumber * _fChargeNumber " << pow(1.+_fRedshift,3.)*dfin* proton_mass_c2 / _fMass / c_squared * _fChargeNumber * _fChargeNumber << std::endl;
	std::cerr << "  (norm_p-lpInt->E_part(index_E) " <<  (norm_p-lpInt->E_part(index_E)) << std::endl;
	std::cerr << " dfin " << dfin << " _fMass " << _fMass << " _fChargeNumber " << _fChargeNumber << " pow(1.+_fRedshift,3.) " << pow(1.+_fRedshift,3.) << std::endl;
	std::cerr << " _fMassNumber " << _fMassNumber << std::endl;
	std::cerr << " _fInitType " << _fInitType << std::endl;
      }
      //std::cout << "AS_TEST PairProduction: " <<  100 * Mpc / _fChargeNumber / _fChargeNumber << std::endl;

    
 
    double lCritStep = 0 ; 

    // If the propagation step is smaller than lCritStep, no secondaries are followed.
    // Allows to avoid secondary propagation when near an observer (step too small)
#define CUT_SECPAIRPROD_SMALLSTEP
#ifdef CUT_SECPAIRPROD_SMALLSTEP
    lCritStep = 0.1*_fpUniverse->IntegratorMinTimeStep() ;
#endif
    if (lpInt->SecPairProdPhotonFlag() && _fTimeStep >= lCritStep &&
	(_fpUniverse->Type() != UNIVERSE_ENV1D  || _fPosition.x() > 0)) {
      double alpha = RandFlat::shoot() ;
      if (alpha <= lpInt->SecPairProdPhotonProba()) {
	TVector3D lPosShower = -_fMomentum.unit() ; // initial shower position in the middle of the step
	lPosShower *= _fTimeStep*c_light/2. ;
	lPosShower += _fPosition ;
	TVector3D lGammaMomentum = _fMomentum ;
	lGammaMomentum.setMag(_fEnergy*loss_pair * inv_c_light) ;
	  if (_fpUniverse->Type() == UNIVERSE_ENV1D) {
	    if (lpInt->KelnerPairProdFlag())
	      //Here a really absurd bug was introduced :) / Upper Line should be commented out !!! JK
	      //_fpList1DPhotons->AddPairProdPhotons(_fPosition.x(), _fEnergy  / ( this->Mass() * c_squared / proton_mass_c2), _fTimeStep * proton_mass_c2 / _fMass / c_squared * _fChargeNumber * _fChargeNumber, _fpUniverse, _fRedshift);
	      _fpList1DPhotons->AddPairProdPhotons(_fPosition.x(), _fEnergy  / ( this->Mass() * c_squared / proton_mass_c2), _fTimeStep * _fChargeNumber * _fChargeNumber, _fpUniverse, _fRedshift);
	    else 
	      _fpList1DPhotons->AddPhotons(_fPosition.x(), _fEnergy, _fEnergy*loss_pair, _fMassNumber) ;
	  } else {
	    if (lpInt->KelnerPairProdFlag())
	      //Second part of the aforementionend bug. 
// 	      fParticleQueue.push_back(new TPhoton(_fpUniverse,_fEnergy,_fTimeStep * proton_mass_c2 / _fMass / c_squared * _fChargeNumber * _fChargeNumber,_fEnergy*loss_pair,lPosShower,
// 						   lGammaMomentum,_fInitPosition,
// 						   _fInitMomentum.getR()*c_light,_fRedshift, _fMassNumber, _fInitType, "PairProd", _fTime));
	      fParticleQueue.push_back(new TPhoton(_fpUniverse,_fEnergy /  ( this->Mass() * c_squared / proton_mass_c2),_fEnergy*loss_pair ,_fTimeStep * _fChargeNumber * _fChargeNumber, lPosShower, lGammaMomentum,_fInitPosition, _fInitMomentum.getR()*c_light,_fRedshift, _fMassNumber, _fInitType, "PairProd", _fTime));
	    else
	      fParticleQueue.push_back(new TPhoton(_fEnergy,_fEnergy*loss_pair,lPosShower,lGammaMomentum,
						   _fInitPosition,_fInitMomentum.getR()*c_light, _fpUniverse, _fMassNumber, _fInitType, "PairProd", _fTime));
	  }

      }//alpha if
    }
  
  _fEnergy = _fMomentum.getR() * c_light ;
    } 
  }
}


double TNucleus::PionProtonMFP( int caseFlag) {
  //Calculates the MFP for protons 
  TInteractionData *lpInt = _fpUniverse->InteractionData() ;
  
  double lRefEnergy = _fEnergy*(1.+_fRedshift) / (double) _fMassNumber; //?
  int i_in;
  double pMFP;
  double scaleFactor=pow(1.+_fRedshift,3);
  //  int lBgFlag;
  //  double lDistance = _fTimeStep * c_light ;
  
  switch (caseFlag){
  case 1:
    //lBgFlag = 1; // CMB
    i_in = min((int)(log10(lRefEnergy/lpInt->E_pionprod(0))/lpInt->dEtabPion()) , 99) ;
    if(i_in == 99){
      pMFP = _fChargeNumber*scaleFactor*lpInt->LossRateProton(i_in); //When the energy is higher than the highest point in the table, the loss rate is set to the loss rate at the highest energy in the table.
      std::cout<< "WARNING: Energy higher than highest energy in pion production table!" << std::endl;
    }
    else
      pMFP = (i_in > 0)?  _fChargeNumber *scaleFactor*
	( (lRefEnergy-lpInt->E_pionprod(i_in)) *(lpInt->LossRateProton(i_in+1)-lpInt->LossRateProton(i_in)) /(lpInt->E_pionprod(i_in+1)-lpInt->E_pionprod(i_in))
	  +lpInt->LossRateProton(i_in)
	  ): 0 ;   
    break;
  case 2:    
    //lBgFlag = 2 ; // IRO
    scaleFactor=_fpUniverse->IRBzEvolutionModel()->GetScalingFactor(_fRedshift);
    i_in = min((int)(log10(lRefEnergy/lpInt->E_IRpionprod(0))/lpInt->dEtabIRPion()) , 149) ;
    if(i_in == 149){
      pMFP = _fChargeNumber*scaleFactor*lpInt->IRLossRateProton(i_in);
      std::cout<< "WARNING: Energy higher than highest energy in pion production table!" << std::endl;
    }
    else
      pMFP = (i_in > 0)?  _fChargeNumber *scaleFactor*
	( (lRefEnergy-lpInt->E_IRpionprod(i_in)) *(lpInt->IRLossRateProton(i_in+1)-lpInt->IRLossRateProton(i_in)) /(lpInt->E_IRpionprod(i_in+1)-lpInt->E_IRpionprod(i_in))
	  +lpInt->IRLossRateProton(i_in)
	  ): 0 ;   
    break;
  case 3:
    if (_fpUniverse->Infrared()->Type() == SHELL)
      pMFP = _fChargeNumber * this->CalcIRLossRate(_fEnergy/GeV, TPhotonSpectrum(_fpUniverse->Infrared()->Spectrum(_fPosition.mag()/Mpc)));
    else 
      pMFP = _fChargeNumber * this->CalcIRLossRate(_fEnergy/GeV, TPhotonSpectrum(_fpUniverse->Infrared()->Spectrum(_fPosition.x()/Mpc,_fPosition.y()/Mpc,_fPosition.z()/Mpc)));
    break;
  default: 
    std::cout<<"Default Case"<< std::endl; 
    throw TCrpErr("TNucleus::PionProtonMFP Unknown photon background for Sophia pion production"); 
  }
  double Nucleicorrection;

  if(_fMassNumber==1)Nucleicorrection=1.;
  else if(_fMassNumber<=8)Nucleicorrection=_fChargeNumber/(0.85*pow(_fChargeNumber,2./3.));
  else if(_fMassNumber>8)Nucleicorrection=1./0.85;

  if(pMFP > 0.)return (Nucleicorrection*inv_c_light / pMFP );
  

   return(1.e99);
}

double TNucleus::PionNeutronMFP( int caseFlag) {
  //Calculates pion production mfp for neutrons

  TInteractionData *lpInt = _fpUniverse->InteractionData() ;

  double lRefEnergy = _fEnergy*(1.+_fRedshift) / (double) _fMassNumber; 
  int i_in;
  double nMFP;
  double scaleFactor=pow(1+_fRedshift,3);
 
  switch (caseFlag){
  case 1:
    //lBgFlag = 1; // CMB
    i_in = min((int)(log10(lRefEnergy/lpInt->E_pionprod(0))/lpInt->dEtabPion()) , 99) ;
    if(i_in == 99){
      nMFP = ( _fMassNumber - _fChargeNumber)*scaleFactor*lpInt->LossRateNeutron(i_in);
      std::cout<< "WARNING: Energy higher than highest energy in pion production table!" << std::endl;
    }
    else
      nMFP = (i_in > 0)?  ( _fMassNumber - _fChargeNumber) *scaleFactor*
	( (lRefEnergy-lpInt->E_pionprod(i_in)) *(lpInt->LossRateNeutron(i_in+1)-lpInt->LossRateNeutron(i_in)) /(lpInt->E_pionprod(i_in+1)-lpInt->E_pionprod(i_in))
	  +lpInt->LossRateNeutron(i_in)
	  ): 0 ;
    break;
  case 2:
    //lBgFlag = 2 ; // IRO
    scaleFactor=_fpUniverse->IRBzEvolutionModel()->GetScalingFactor(_fRedshift);//IRO scaling
    i_in = min((int)(log10(lRefEnergy/lpInt->E_IRpionprod(0))/lpInt->dEtabIRPion()) , 149) ;
    if(i_in == 149){
      nMFP = ( _fMassNumber - _fChargeNumber)*scaleFactor*lpInt->IRLossRateNeutron(i_in);
      std::cout<< "WARNING: Energy higher than highest energy in pion production table!" << std::endl;
    }
    else
      nMFP = (i_in > 0)?  ( _fMassNumber - _fChargeNumber) *scaleFactor*
	( (lRefEnergy-lpInt->E_IRpionprod(i_in)) *(lpInt->IRLossRateNeutron(i_in+1)-lpInt->IRLossRateNeutron(i_in)) /(lpInt->E_IRpionprod(i_in+1)-lpInt->E_IRpionprod(i_in))
	  +lpInt->IRLossRateNeutron(i_in)
	  ): 0 ;   
    break;
  case 3:
    if (_fpUniverse->Infrared()->Type() == SHELL)
      nMFP = ( _fMassNumber - _fChargeNumber) * this->CalcIRLossRate(_fEnergy/GeV, TPhotonSpectrum(_fpUniverse->Infrared()->Spectrum(_fPosition.mag()/Mpc)));
    else 
      nMFP = ( _fMassNumber - _fChargeNumber) * this->CalcIRLossRate(_fEnergy/GeV, TPhotonSpectrum(_fpUniverse->Infrared()->Spectrum(_fPosition.x()/Mpc,_fPosition.y()/Mpc,_fPosition.z()/Mpc)));
    
    break;
  default: 
    throw TCrpErr("TNucleus::PionNeutronMFP Unknown photon background for Sophia pion production"); 
  }
  double Nucleicorrection;

  if(_fMassNumber==1)Nucleicorrection=1.;
  else if(_fMassNumber<=8)Nucleicorrection=(_fMassNumber - _fChargeNumber)/(0.85*pow(_fMassNumber - _fChargeNumber,2./3.));
  else if(_fMassNumber>8)Nucleicorrection=1./0.85;

  if(nMFP >0.)return (Nucleicorrection*inv_c_light /  nMFP);
  return(1.e99);
}


void TNucleus::SophiaPionProd(const int aBgFlag, 
			      int aNucleonChargeNumber,
			      double EnergyPerNucleon) {

  // p + g -> p + Pi0 
  // p + g -> n + Pi+
  // n + g -> n + Pi0
  // n + g -> p + Pi-

  //std::cout<<"TEST : Start SophiaPionProd()"<<std::endl;

  TInteractionData *lpInt = _fpUniverse->InteractionData() ;
  
  int lNature = 0;

  static int testflag=0; 
  if(testflag==0){
    testflag=1;
    std::cout<< "\t->PionProduction called" << std::endl; 
    //std::cout<< "Achtung modified!!! " << std::endl; 
  }


  //no anti-nucleons supported right now! To make surce
  if(_fCharge < 0. || aNucleonChargeNumber < 0.) throw TCrpErr("Error: Nucleus with negative Charge! SophiaPionProd() doesn't support that!"); // TODO : precuation : check if this is needed.
  

  if (aNucleonChargeNumber == 0) lNature = 1 ;       //incident neutron
  double lEinGeV = EnergyPerNucleon * invGeV ;
  //std::cout<<"TEST EnergyPerNucleon = "<<EnergyPerNucleon<<std::endl;
  double lMomentaList[5][2000] ; 
  int lParticleList[2000] ;
  int lNbParticle, lPartType ;
  int lCheckFlag= 0 ;                                     // to check if only one p or n is out!
  int lBgFlag = aBgFlag ;
  double lRedshift = _fRedshift ;
  double lZmax = _fpUniverse->Infrared()->Zmax() ;  
  /*
  if (lBgFlag != 3)
    sophiaevent(lNature, lEinGeV, lMomentaList, lParticleList, lNbParticle, 
		lRedshift, lBgFlag, lZmax) ;
  else {
  */

  //  cout << "In SophiaPionProd with lBgFlag = " << lBgFlag << endl;

  int lIRNE = 2;
  double lDummy[2];
  if (lBgFlag < 3) {
#ifdef DEBUG_OUTPUT 
    std::cout<< "Vor sophia call" << std::endl;
#endif 
    sophiaevent_(lNature, lEinGeV, lMomentaList, lParticleList, lNbParticle, 
		      lRedshift, lBgFlag, lZmax, lIRNE, lDummy, lDummy);    
#ifdef DEBUG_OUTPUT 
    std::cout<< "Nach sophia call" << std::endl;
#endif 
  }
  else {
    lIRNE = _fpUniverse->Infrared()->NE();
    double lIRenergy[lIRNE];
    double lIRspectrum[lIRNE];
    TPhotonSpectrum phsp;
    if (_fpUniverse->Infrared()->Type() == SHELL)
      phsp = _fpUniverse->Infrared()->Spectrum(_fPosition.mag()/Mpc);
    else
      phsp = _fpUniverse->Infrared()->Spectrum(_fPosition.x()/Mpc,_fPosition.y()/Mpc,_fPosition.z()/Mpc);

    for (int i = 0; i< lIRNE; i++) {
      lIRenergy[i] = phsp.GetE()[i];
      lIRspectrum[i] = phsp.Spectrum()[i];
    }
    
    sophiaevent_(lNature, lEinGeV, lMomentaList, lParticleList, lNbParticle, 
		      lRedshift, lBgFlag, lZmax, lIRNE, lIRenergy, lIRspectrum);    
  }
  // important : the incoming momentum is along z for sophia.
  // we do only use the outcoming energy, ie. lMomentaList[3][]
  for (int i=0; i<lNbParticle; i++) {
    lPartType = lParticleList[i] ;
    TVector3D lMomentum = _fMomentum;
    lMomentum.setMag(lMomentaList[3][i] * GeV * inv_c_light);
    //std::cout<<"TEST Energy of the particle"<<lMomentaList[3][i] * GeV<<std::endl;
    if (lPartType == 13 || lPartType == 14) {                  
      //proton or neutron has interacted -> change mass and charge number + add the then free
      //nucleon to the particle que.
      if (lCheckFlag == 0) {
	lCheckFlag += 1 ;
	//the nucleon which has interacted leaves the nucleus -> modify charge- and massnumber.
	if (lPartType == 13){                                   //proton after pion prod.
	  if(aNucleonChargeNumber == 1){                        //p->p
	    //	    std::cout<<"p->p"<<std::endl;
	    if(_fMassNumber >= 2){                              
	      _fChargeNumber-- ;
	      _fMassNumber-- ;
	      _fMass   = this->Mass();
	      _fCharge = eplus *_fChargeNumber ;   
	      _fEnergy -= EnergyPerNucleon;
	      _fMomentum.setMag(_fEnergy * inv_c_light) ; 
	      fParticleQueue.push_back(new TNucleus(_fpUniverse, _fPosition, lMomentum, 
						    _fInitPosition, _fInitMomentum, 
						    _fTime, 1, 1, _fpList1DPhotons, _fInitType
						    )) ;
	    }else{	     
	      _fEnergy = lMomentaList[3][i] * GeV;
	      _fMomentum.setMag(_fEnergy * inv_c_light) ; 
	    }
	  }else if(aNucleonChargeNumber == 0){                  //n->p     
	    //std::cout<<"Test neutron reacts n->p"<<std::endl;
	    if(_fMassNumber >= 2){
	      _fMassNumber--;
	      _fMass   = this->Mass() ;
	      _fEnergy -= EnergyPerNucleon;
	      _fMomentum.setMag(_fEnergy * inv_c_light) ; 
	      fParticleQueue.push_back(new TNucleus(_fpUniverse, _fPosition, lMomentum, 
						    _fInitPosition, _fInitMomentum, 
						    _fTime, 1, 1, _fpList1DPhotons, _fInitType
						    )) ;  
	    }else{	     
	      _fChargeNumber=1 ;
	      _fMass   = this->Mass() ;
	      _fCharge = eplus *_fChargeNumber ;
	      _fEnergy = lMomentaList[3][i] * GeV;
	      _fMomentum.setMag(_fEnergy * inv_c_light) ; 
	    }
	    
	  }else throw TCrpErr("Error in SophiaPionProd(): Unknown input particle.");	  
	}
	else if (lPartType == 14){                              //neutron after pion prod.  
	  //std::cout<<"TEST: neutron created in PiP!"<<std::endl;
	  if(aNucleonChargeNumber == 1){                        //p->n
	    //	    std::cout<<"p->n"<<std::endl;
	    if(_fMassNumber >= 2){
	      _fChargeNumber-- ;
	      _fMassNumber-- ;
	      _fMass   = this->Mass() ;
	      _fCharge = eplus *_fChargeNumber ;
	      _fEnergy -= EnergyPerNucleon;
	      _fMomentum.setMag(_fEnergy * inv_c_light) ; 
	      fParticleQueue.push_back(new TNucleus(_fpUniverse, _fPosition, lMomentum, 
						    _fInitPosition, _fInitMomentum, 
						    _fTime, 1, 0, _fpList1DPhotons, _fInitType
						    )) ;


	    }else{
	      _fChargeNumber=0;
	      _fMass   = this->Mass() ;
	      _fCharge = 0. ;
	      _fEnergy = lMomentaList[3][i] * GeV;
	      _fMomentum.setMag(_fEnergy * inv_c_light) ; 
	    }
	  }else if(aNucleonChargeNumber == 0){                 //n->n     
	    //	    std::cout<<"n->n"<<std::endl;
	    if(_fMassNumber >= 2){
	      _fMassNumber--;
	      _fMass   = this->Mass() ;
	      _fEnergy -= EnergyPerNucleon;
	      _fMomentum.setMag(_fEnergy * inv_c_light) ; 
	      fParticleQueue.push_back(new TNucleus(_fpUniverse, _fPosition, lMomentum, 
						    _fInitPosition, _fInitMomentum, 
						    _fTime, 1, 0, _fpList1DPhotons , _fInitType
						    )) ;
	    }else{
	      _fEnergy = lMomentaList[3][i] * GeV;
	      _fMomentum.setMag(_fEnergy * inv_c_light) ; 
	    }
	  }else throw TCrpErr("Error in SophiaPionProd(): Unknown input particle.");
	}
      } else { // multinucleus production --> need to add a nucleus in the queue
	double lCharge = 0 ;
	if (lPartType == 13) lCharge = 1 ;
	if (lPartType == 14) lCharge = 0 ;
	if(lPartType == 13) pprod +=1;
	if(lPartType == 14) nprod +=1;
	TVector3D lMomentum = _fMomentum ;
	lMomentum.setMag( lMomentaList[3][i] * GeV * inv_c_light ) ;
	fParticleQueue.push_back(new TNucleus(_fpUniverse, _fPosition, lMomentum, _fInitPosition, 
					      _fInitMomentum, _fTime, 1, lCharge, _fpList1DPhotons, _fInitType
						    )) ;
      }
    } else if (lPartType == 1 || lPartType == 2 || lPartType == 3) {
      // New electromagnetic cascade.
      double lCharge = 0;
      if (lPartType == 2) lCharge = eplus ;
      if (lPartType == 3) lCharge = electron_charge ;
      TVector3D lGammaMomentum = _fMomentum ;
      lGammaMomentum.setMag(lMomentaList[3][i] * GeV * inv_c_light) ;
      if (lpInt->SecPhotonFlag()) 
	if (_fpUniverse->Type() == UNIVERSE_ENV1D) {
	  PARTICLE lPart ;
	  if (lCharge == 0) {
	    lPart = PHOTON ;
	  }else if (lCharge == eplus) {
	    lPart = POSITRON ;
	  } else if (lCharge == electron_charge) {
	    lPart = ELECTRON ;
	  } else throw TCrpErr("Charge error!") ;
	  _fpList1DPhotons->AddPhotons(_fPosition.x(), lGammaMomentum.mag()*c_light, lPart) ;
	} else { // Universe Largescale
	  if(lEinGeV*GeV < lGammaMomentum.getR()*c_light) {
	    std::cout << " (lEinGeV < lGammaMomentum.getR()*c_light)\n _fEnergy "<< lEinGeV << std::endl; 
	    std::cout << " A " << _fMassNumber << "  Z " << _fChargeNumber << std::endl; 
	    std::cout << " _fEnergy " << _fEnergy << " lGammaMomentum.mag()*c_light " << lGammaMomentum.mag()*c_light << std::endl;
	    std::cout << " _fInitMomentum.getR()*c_light " << _fInitMomentum.getR()*c_light << std::endl; 
	    exit(-1);
	  }
	  string lPart ;
	  if (lCharge == 0) {
	    lPart = "pi_gamma";
	  }else if (lCharge == eplus) {
	    lPart = "pi_e+";
	  } else if (lCharge == electron_charge) {
	    lPart = "pi_e-" ;
	  } else throw TCrpErr("Charge error!") ;
	  fParticleQueue.push_back(new TPhoton(_fPosition, lGammaMomentum, lCharge, _fInitPosition, 
					       _fInitMomentum.getR()*c_light, _fpUniverse, _fInitType, lPart, _fTime));
	}
    } else if (lPartType >= 15 && lPartType <= 18) {
      // New neutrino!
      int lFlavor = 12 ; 
      if (lPartType == 15) lFlavor = 12 ; // nu_e
      if (lPartType == 16) lFlavor = -12 ; // nu_e~
      if (lPartType == 17) lFlavor = 14 ; // nu_mu
      if (lPartType == 18) lFlavor = -14 ; // nu_mu~
      TVector3D lNuMomentum = _fMomentum ;
      lNuMomentum.setMag(lMomentaList[3][i] * GeV * inv_c_light) ;
      if (lpInt->SecNuFlag()) 
	fParticleQueue.push_back(new TNeutrino(lFlavor,_fPosition, lNuMomentum, _fInitPosition, _fInitMomentum.mag()*c_light, _fpUniverse, _fInitType, _fTime)) ;
      
    } else if (lPartType == -13 || lPartType == -14) {
      // Antiproton/neutron --> forget it! (very small flux)
      bbarprod+=1;
    } else throw TCrpErr("Unknown particle generated by Sophia.") ;

  }
  //  if (lCheckFlag == 0) throw TCrpErr("No nucleus generated by Sophia.") ;
  if (lCheckFlag == 0) std::cout << " WARNING: No nuclues generated by Sophia, new interaction dialed." << std::endl;

}


inline double TNucleus::Mass(){
  //Just a wrapper for TNucleusDB::GetMass()
  return TNucleusDB::GetInstance()->GetMass(_fChargeNumber * 1000 + _fMassNumber);
}  

void TNucleus::BetaMinusDecay(){
  _fChargeNumber+=1;

  double E, lQ, lGamma;
  lGamma = _fEnergy / _fMass / c_squared;
  lQ = _fMass;
  if (lQ > this->Mass()){ 
    lQ -= this->Mass(); 
  } else {// should never be triggered!
    std::cout<< "in exception path"<< std::endl;
    std::cout<<  TNucleusDB::GetInstance()->GetMass(_fChargeNumber * 1000 + _fMassNumber)<< std::endl;
    std::cout<< "Fe Mass: " << TNucleusDB::GetInstance()->GetMass(26056)/c_squared << std::endl;
    std::cout<< "Old Mass*c^2 : " << TNucleusDB::GetInstance()->GetMass(_fChargeNumber-1 * 1000 + _fMassNumber)*c_squared << std::endl;
    std::cout<< "New Mass*c^2 : " << TNucleusDB::GetInstance()->GetMass(_fChargeNumber * 1000 + _fMassNumber)*c_squared << std::endl;
    _fEnergy = lQ; 
    std::cout<< " WARNING: Decay : lQ < _fMass " << std::endl;
    std::cout<< "A " << _fMassNumber << "  Z  " << _fChargeNumber << "  Decay B-" << std::endl; 
    throw TCrpErr("TNucleus::BetaMinusDecay Decay not allowed!!");
  }//Checks if decay is energetically allowed. 
  lQ*=c_squared;
  // std::cout<< " lQ " << lQ << std::endl; 
  _fEnergy *=this->Mass() / _fMass; // nucleus rests in CMS
  _fMass=this->Mass(); 
  _fCharge = eplus * _fChargeNumber ;
  _fMomentum.setMag(_fEnergy / c_light) ;
  // E = _fEnergy * invEeV ; // f77 units 


  TInteractionData *lpInt = _fpUniverse->InteractionData() ;
  //#define betaSpec
#ifdef betaSpec
  fstream betaOut ("betaSpec.txt", fstream::out);
  if( betaOut.fail()) throw TCrpErr("betaOut Fail!");
  betaOut << "#Beta spectrum of A=" << _fMassNumber << " Z=" << _fChargeNumber << std::endl;
  betaOut << "#lQ=" << lQ << " electron_mass_c2="<< electron_mass_c2 << std::endl;
  betaOut << "#randumNumber  lEnergy_e_rest " << std::endl;
  for(int i=0; i<10000; i++){
    double outlEnergy_e_rest = lpInt->RandDistriNeutronDecay()->shoot();
    betaOut << outlEnergy_e_rest << "\t" << electron_mass_c2 + outlEnergy_e_rest*(lQ-electron_mass_c2) << std::endl;
  }
  betaOut.close();
  throw TCrpErr("STOP! betaSpec.txt created.");
#endif

    // Secondary electron and antineutrino

    if ((lpInt->SecPhotonFlag() || lpInt->SecNuFlag()) &&
	(_fpUniverse->Type() != UNIVERSE_ENV1D  || _fPosition.x() > 0)) {
      double lEnergy_e_rest = lpInt->RandDistriNeutronDecay()->shoot() ; // Gute naeherung? Nein, S.202 im Basdevant
      lEnergy_e_rest = electron_mass_c2 + lEnergy_e_rest*(lQ-electron_mass_c2) ; 
      if (lpInt->SecPhotonFlag()) {
	double lP_e_rest = sqrt(max(lEnergy_e_rest*lEnergy_e_rest-electron_mass_c2*electron_mass_c2,0.)) ;
	double lCosTheta_e_rest = 2*RandFlat::shoot()-1 ; 
	double lEnergy_e = lGamma * (lEnergy_e_rest + lP_e_rest*lCosTheta_e_rest) ; // + because cos theta = - cos (theta + pi)
	TVector3D lGammaMomentum = _fMomentum ;
	lGammaMomentum.setMag(lEnergy_e * inv_c_light ) ;
	if (_fpUniverse->Type() == UNIVERSE_ENV1D) {
	  PARTICLE lPart = ELECTRON ;
	  _fpList1DPhotons->AddPhotons(_fPosition.x(), lGammaMomentum.mag()*c_light, lPart) ;
	  if( lGammaMomentum.mag()*c_light > _fEnergy){
	    std::cout << "Beta- Decay,  lGammaMomentum.mag()*c_light > _fEnergy" << std::endl; 
	    std::cout << " A  " << _fMassNumber << " Z " << _fChargeNumber << std::endl; 
	    std::cout << " _fEnergy " << _fEnergy << "  E = _fEnergy * invEeV " << E << std::endl; 
	    std::cout << " lGammaMomentum.mag()*c_light " << lGammaMomentum.mag()*c_light <<std::endl;
	    std::cout << " lQ " << lQ << " old lQ " << neutron_mass_c2 - proton_mass_c2 << std::endl;
	    std::cout << " lP_e_rest " << lP_e_rest <<std::endl;
	    exit(-1); 
	  }
	} else{
	  char strg[11];
	  sprintf(strg, "A:%iZ:%ib-\0", _fMassNumber, _fChargeNumber - 1);
	  fParticleQueue.push_back(new TPhoton(_fPosition, lGammaMomentum, -1 * electron_charge,
					       _fInitPosition, _fInitMomentum.getR()*c_light, _fpUniverse, _fInitType, strg, _fTime));
	}
      }// SecPhotonFlag
      if (lpInt->SecNuFlag()) {
	double lEnergy_nu_rest = lQ - lEnergy_e_rest ;
	double lCosTheta_nu_rest = 2*RandFlat::shoot()-1 ;
	double lEnergy_nu = lGamma * lEnergy_nu_rest * (1. + lCosTheta_nu_rest) ; // + because cos theta = - cos (theta + pi)
	TVector3D lNuMomentum = _fMomentum ;
	lNuMomentum.setMag(lEnergy_nu * inv_c_light ) ;
	fParticleQueue.push_back(new TNeutrino(-12, _fPosition, lNuMomentum, _fInitPosition, _fInitMomentum.mag()*c_light, _fpUniverse, _fInitType, _fTime)) ; // nu_e~ flavor code
      }//SecNuFlag
    }//Long if
}

void TNucleus::BetaPlusDecay(){
  _fChargeNumber+=-1;
  
  double E, lQ, lGamma;
  lGamma = _fEnergy / _fMass / c_squared;
  lQ = _fMass;
  if (lQ > this->Mass()){ 
    lQ -= this->Mass(); 
  } else {// should never be triggered!
    _fEnergy = lQ; 
    
    std::cout<< " WARNING: Decay : lQ < _fMass " << std::endl;
    std::cout<< "Fe Mass: " << TNucleusDB::GetInstance()->GetMass(26056)/c_squared << std::endl;
    std::cout<< "Old Mass/c^2 : " << TNucleusDB::GetInstance()->GetMass(_fChargeNumber+1 * 1000 + _fMassNumber) << std::endl;
    std::cout<< "New Mass/c^2 : " << TNucleusDB::GetInstance()->GetMass(_fChargeNumber * 1000 + _fMassNumber) << std::endl;
    std::cout<< "A " << _fMassNumber << "  Z  " << _fChargeNumber << "  Decay B+ " << std::endl; 
    throw TCrpErr("TNucleus::GetNuclearDecayProduct Decay not allowed!!");
  }//Checks if decay is energetically allowed. 
  lQ*=c_squared;
  // std::cout<< " lQ " << lQ << std::endl; 
  _fEnergy *=this->Mass() / _fMass; // nucleus rests in CMS
   
  if ( lQ < electron_mass_c2) {
    std::cout << " only electron capture allowed, A: " << _fMassNumber << " Z=" << _fChargeNumber <<  " DecayType B+" << std::endl; 
    throw TCrpErr("TNucleus::BetaPlusDecay Decay not allowed!!");
  }
  _fMass=this->Mass(); 
  _fCharge = eplus * _fChargeNumber ;
  _fMomentum.setMag(_fEnergy / c_light) ;
  // E = _fEnergy * invEeV ; // f77 units 
    // Secondary positron and neutrino
    TInteractionData *lpInt = _fpUniverse->InteractionData() ;
    if ((lpInt->SecPhotonFlag() || lpInt->SecNuFlag()) &&
	(_fpUniverse->Type() != UNIVERSE_ENV1D  || _fPosition.x() > 0)) {
      double lEnergy_e_rest = lpInt->RandDistriNeutronDecay()->shoot() ; // Gute naeherung? Nein, S.202 im Basdevant
      lEnergy_e_rest = electron_mass_c2 + lEnergy_e_rest*(lQ-electron_mass_c2) ; 
      if (lpInt->SecPhotonFlag()) {
	double lP_e_rest = sqrt(max(lEnergy_e_rest*lEnergy_e_rest-electron_mass_c2*electron_mass_c2,0.)) ;
	double lCosTheta_e_rest = 2*RandFlat::shoot()-1 ; 
	double lEnergy_e = lGamma * (lEnergy_e_rest + lP_e_rest*lCosTheta_e_rest) ; //SOLLTE DAS MINUS SEIN???? 
	TVector3D lGammaMomentum = _fMomentum ;
	lGammaMomentum.setMag(lEnergy_e * inv_c_light ) ;
	//std::cout<<"\n\n\t\t\t TEST2m!  \n\n\n"<<std::endl;
	if (_fpUniverse->Type() == UNIVERSE_ENV1D) {
	  PARTICLE lPart = POSITRON ;
	  if( lGammaMomentum.mag()*c_light > _fEnergy){
	    std::cout << "Beta+ Decay,  lGammaMomentum.mag()*c_light > _fEnergy" << std::endl; 
	    std::cout << " A  " << _fMassNumber << " Z " << _fChargeNumber << std::endl; 
	    std::cout << " _fEnergy " << _fEnergy << "  E = _fEnergy * invEeV " << E << std::endl; 
	    std::cout << " lGammaMomentum.mag()*c_light " << lGammaMomentum.mag()*c_light <<std::endl;
	    exit(-1); 
	  }
	  _fpList1DPhotons->AddPhotons(_fPosition.x(), lGammaMomentum.mag()*c_light, lPart) ;
	} else{
	  char strg[11];
	  //std::cout<<"\n\n\t\t\t TEST2n!  \n\n\n"<<std::endl;
	  sprintf(strg, "A:%iZ:%ib+\0", _fMassNumber, _fChargeNumber + 1);
	  //std::cout<<"\n\n\t\t\t TEST2o!  \n\n\n"<<std::endl;
	  fParticleQueue.push_back(new TPhoton(_fPosition, lGammaMomentum, electron_charge,
					       _fInitPosition, _fInitMomentum.getR()*c_light, _fpUniverse, _fInitType, strg, _fTime));
	}
      }// SecPhotonFlag
      if (lpInt->SecNuFlag()) {
	double lEnergy_nu_rest = lQ - lEnergy_e_rest ;
	double lCosTheta_nu_rest = 2*RandFlat::shoot()-1 ;
	double lEnergy_nu = lGamma * lEnergy_nu_rest * (1. + lCosTheta_nu_rest) ;
	TVector3D lNuMomentum = _fMomentum ;
	lNuMomentum.setMag(lEnergy_nu * inv_c_light ) ;
	fParticleQueue.push_back(new TNeutrino(12, _fPosition, lNuMomentum, _fInitPosition, _fInitMomentum.mag()*c_light, _fpUniverse, _fInitType, _fTime)) ; // nu_e~ flavor code
      }//SecNuFlag
    }//Long if
}

void TNucleus::AlphaDecay(){
  //Change of the nucleus
  TVector3D pMomentum= _fMomentum;
  double lGamma;
  _fChargeNumber+=-2;
  _fMassNumber+=-4;
  lGamma = _fEnergy / _fMass / c_squared;
  _fEnergy *=this->Mass() / _fMass; // nucleus rests in CMS
  _fMass=this->Mass();
  _fCharge = eplus * _fChargeNumber ;
  _fMomentum.setMag(_fEnergy / c_light) ;
  
  //Creation of the alpha particle
  pMomentum.setMag(lGamma * TNucleusDB::GetInstance()->GetMass(2004)*c_light) ;
  fParticleQueue.push_back(new TNucleus(_fpUniverse, _fPosition, pMomentum, _fInitPosition, _fInitMomentum, _fTime, 4, 2,_fpList1DPhotons, _fInitType));
}


void TNucleus::ProtonDripping(){
  TVector3D pMomentum= _fMomentum;
  double lGamma;
  _fChargeNumber+=-1;
  _fMassNumber+=-1;

  lGamma = _fEnergy / _fMass / c_squared;
  _fEnergy *=this->Mass() / _fMass; // nucleus rests in CMS
  _fMass=this->Mass();
  _fCharge = eplus * _fChargeNumber ;
  _fMomentum.setMag(_fEnergy / c_light) ;
  
  //Creation of the secondary proton
    pMomentum.setMag(lGamma * proton_mass_c2 * inv_c_light) ;
    fParticleQueue.push_back(new TNucleus(_fpUniverse, _fPosition, pMomentum, _fInitPosition, _fInitMomentum, _fTime, 1, 1,_fpList1DPhotons, _fInitType )) ;
}

void TNucleus::NeutronDripping(){
  TVector3D pMomentum= _fMomentum;
  double lGamma;
  _fMassNumber+=-1;

  lGamma = _fEnergy / _fMass / c_squared;
  _fEnergy *=this->Mass() / _fMass; // nucleus rests in CMS
  _fMass=this->Mass();
  _fCharge = eplus * _fChargeNumber ;
  _fMomentum.setMag(_fEnergy / c_light) ;

    pMomentum.setMag(lGamma * neutron_mass_c2 * inv_c_light);
    //  std::cout<< "#Particels" <<fParticleQueue.size()<<std::endl;
    fParticleQueue.push_back(new TNucleus(_fpUniverse, _fPosition, pMomentum, _fInitPosition, _fInitMomentum, _fTime, 1, 0,_fpList1DPhotons, _fInitType )) ;
}


void TNucleus::GetNuclearDecayProducts(){

  //List of Decaytypes: 
  //     "B-" :  "B"
  //     "B+" :  "E" (from electron capture)
  //     "alpha" : "A"
  //     "Proton dripping" : "P"
  //     "Neutron dripping" : "N" 
  //None 0

  const char* decaytype=TNucleusDB::GetInstance()->GetDecayMode(_fChargeNumber * 1000 + _fMassNumber, RandFlat::shoot()*100 ).c_str();

  static int testflag=0; 
  if(testflag==0){
    testflag=1;
    std::cout<< "\t->Decay called" << std::endl; 
  }

  while(*decaytype !='\0'){// NULL terminated c string
    switch(*decaytype){
    case 'B':
      BetaMinusDecay();
      //  std::cout << "BetaDecay" << std::endl;
      break;
    case 'E':
      BetaPlusDecay();
      //  std::cout << "Beta+Decay" << std::endl;
      break;
    case 'P':
      ProtonDripping();
      //  std::cout << "ProtonDripping" << std::endl;
      break;
    case 'N': 
      NeutronDripping();
      //  std::cout << "NeutronDripping" << std::endl;
      break;
    case 'A':
      AlphaDecay();
      //  std::cout << "AlphaDecay" << std::endl;
      break;
    default: 
      std::cout<< "Decay: " << *decaytype << std::endl;
      throw TCrpErr("TNucleus::GetNuclearDecayProduct Unknown Decay.");
    };
    decaytype++;
  }
}


double TNucleus::DecayGammaTau() {
  double tau;
  
  /*Caching */
  static int lastChargeNumber=-1;
  static int lastMassNumber=-1;
  static double lasttau=-1.;
  
  if(lastChargeNumber==_fChargeNumber && lastMassNumber == _fMassNumber) {
    
    return lasttau*_fEnergy / (_fMass) / c_squared;
  }

  lastChargeNumber=_fChargeNumber;
  lastMassNumber=_fMassNumber;

  tau= TNucleusDB::GetInstance()->GetTau(_fChargeNumber * 1000 + _fMassNumber);
  if(tau ==0) tau=1.e99;//Nucleus is stable.
  lasttau=tau;
  tau *=  _fEnergy / (_fMass) / c_squared;
  return tau;
}


void TNucleus::TALYSNucleiPhotoDisintegration(TBasicParam* aBasic){
  
 static int testflag=0; 
 if(testflag==0){
  testflag=1;
  std::cout<< "\t->PhotoDisintegration called" << std::endl; 
 }

  TInteractionData *lpInt = _fpUniverse->InteractionData();

  //Note, effect of redshift on mean free path is l(E, z)^-1 = (1 + z)^3 * l[(1 + z)E, z = 0]
  //see CRPropa paper Armengaud et al eq. (13).
  double lRedshiftedEnergy = _fEnergy*(1.+_fRedshift);

  //Get exclusive channel 
  unsigned long int InitNucleusId = _fChargeNumber*1000+_fMassNumber;

  unsigned long int TheExclusiveChannel = 
    lpInt->GetExclusivePDChannelDirectly(InitNucleusId,
				       _fPosition,
				       lRedshiftedEnergy/1000.,  //GeV
				       _fMass*c_squared/1000.,   //GeV
				       _fRedshift,
				       _fTimeStep);
   double lEnergyPerNucleon = _fEnergy/((double) _fMassNumber);
  
  //change the properties of the photodisintegrated "mother" nucleus
  std::vector<int> DMassAndDCharge;
  DMassAndDCharge=(lpInt->GetPointerToPhotoDisTables())[0]->GetDMassAndDChargeExclId(TheExclusiveChannel); 
  
  int DeltaMassNumber   = DMassAndDCharge[0]; 
  int DeltaChargeNumber = DMassAndDCharge[1];

  _fMassNumber   -= DeltaMassNumber;
  _fChargeNumber -= DeltaChargeNumber;
  _fMass          = this->Mass();
  _fCharge        = eplus * (double) _fChargeNumber;
  _fEnergy        = lEnergyPerNucleon * (double) _fMassNumber;
  
  _fMomentum.setMag(_fEnergy * inv_c_light);
  
  //Put the created particles on the stack
  int NNeutron         = (int) (TheExclusiveChannel/1.e5);
  TheExclusiveChannel -= (long unsigned int) (NNeutron*1.e5);
  for(int i=0; i<NNeutron; i++){
    if(lEnergyPerNucleon < aBasic->Emin()) continue;
    TVector3D lMomentum;
    lMomentum=_fMomentum;
    lMomentum.setMag(lEnergyPerNucleon * inv_c_light);
    fParticleQueue.push_back(new TNucleus(_fpUniverse, 
					  _fPosition, 
					  lMomentum, 
					  _fInitPosition, 
					  _fInitMomentum, 
					  _fTime, 
					  1, 
					  0, 
					  _fpList1DPhotons, _fInitType
						    ));
  }
  
  int NProton         = (int) (TheExclusiveChannel/1.e4);
  TheExclusiveChannel-= (long unsigned int) (NProton*1.e4);
  for(int i=0; i<NProton; i++){
    if(lEnergyPerNucleon < aBasic->Emin()) continue;
    TVector3D lMomentum;
    lMomentum=_fMomentum;
    lMomentum.setMag(lEnergyPerNucleon * inv_c_light);
    fParticleQueue.push_back(new TNucleus(_fpUniverse, 
					  _fPosition, 
					  lMomentum, 
					  _fInitPosition, 
					  _fInitMomentum, 
					  _fTime, 
					  1, 
					  1, 
					  _fpList1DPhotons, _fInitType
						    ));
  }
  
  int NDeuterium      = (int) (TheExclusiveChannel/1.e3);
  TheExclusiveChannel-= (long unsigned int) (NDeuterium*1.e3);
  for(int i=0; i<NDeuterium; i++){
    if(lEnergyPerNucleon * 2. < aBasic->Emin()) continue;
    TVector3D lMomentum;
    lMomentum=_fMomentum;
    lMomentum.setMag(lEnergyPerNucleon * 2. * inv_c_light);
    fParticleQueue.push_back(new TNucleus(_fpUniverse, 
					  _fPosition, 
					  lMomentum, 
					  _fInitPosition, 
					  _fInitMomentum, 
					  _fTime, 
					  2, 
					  1, 
					  _fpList1DPhotons, _fInitType
					  ));
  }
  
  int NTritium        = (int) (TheExclusiveChannel/1.e2);
  TheExclusiveChannel-= (long unsigned int) (NTritium*1.e2);
  for(int i=0; i<NTritium; i++){
    if(lEnergyPerNucleon * 3. < aBasic->Emin()) continue;
    TVector3D lMomentum;
    lMomentum=_fMomentum;
    lMomentum.setMag(lEnergyPerNucleon * 3. * inv_c_light);
    fParticleQueue.push_back(new TNucleus(_fpUniverse, 
					  _fPosition, 
					  lMomentum, 
					  _fInitPosition, 
					  _fInitMomentum, 
					  _fTime, 
					  3, 
					  1, 
					  _fpList1DPhotons, _fInitType
					  ));
  }
  
  int NHelium         = (int) (TheExclusiveChannel/1.e1);
  TheExclusiveChannel-= (long unsigned int) (NHelium*1.e1);
  for(int i=0; i<NHelium; i++){
    if(lEnergyPerNucleon * 3. < aBasic->Emin()) continue;
    TVector3D lMomentum;
    lMomentum=_fMomentum;
    lMomentum.setMag(lEnergyPerNucleon * 3. * inv_c_light);
    fParticleQueue.push_back(new TNucleus(_fpUniverse, 
					  _fPosition, 
					  lMomentum, 
					  _fInitPosition, 
					  _fInitMomentum, 
					  _fTime, 
					  3, 
					  2, 
					  _fpList1DPhotons, _fInitType
						    ));
  }
  
  int NAlpha          = (int) (TheExclusiveChannel);
  for(int i=0; i<NAlpha; i++){
    if(lEnergyPerNucleon * 4. < aBasic->Emin()) continue;
    TVector3D lMomentum;
    lMomentum=_fMomentum;
    lMomentum.setMag(lEnergyPerNucleon * 4. * inv_c_light);
    fParticleQueue.push_back(new TNucleus(_fpUniverse, 
					  _fPosition, 
					  lMomentum, 
					  _fInitPosition, 
					  _fInitMomentum, 
					  _fTime, 
					  4, 
					  2, 
					  _fpList1DPhotons, _fInitType
					  ));
  }
}

// Start variable backgrounds part.

void TNucleus::CreatePionTable(){
  TUniformInfrared iro;
  IR iro2;
  std::cout<< "1"  << std::endl;
  ofstream outFile( "pionTable.dat"); //TestCODE
    if( outFile.fail()) throw TCrpErr("CreatePionTable: No outfile");

  double outE;
  std::cout << "Vor vector:" << std::endl;
  vector<double>XData, YData;
  for(double eps=1e-3; eps <= 100.5; eps*=pow(10., 1/3.)) {
    std::cout << "eps " << eps << std::endl;
    XData.push_back(eps ); //Energy in eV 
    YData.push_back(iro.GetPhotonDensity(1., 1., 1., 0, eps ) * pow(eps, 2));
  }
  ofstream outSpec("OutSpec.dat");
  for(int i =0; i < XData.size(); i++) outSpec << XData[i] << "  " << YData[i] << std::endl;
  outSpec.close();
  std::cout << "vor constructor" << std::endl;
  TPhotonSpectrum Spec(XData, YData);
  std::cout << "nach constructor" << std::endl;
  for(outE = 1.e10 * MeV; outE < 1.e18 * MeV; outE *= pow(10., 1./3.)){
    std::cout << outE/ EeV << ':' <<  CalcIRLossRate(outE / GeV ,  Spec) * Mpc << std::endl;
    outFile << outE/ EeV << ' ' <<  CalcIRLossRate(outE / GeV ,  Spec) * Mpc << std::endl;
    // outFile << outE/ EeV << CalcIRLossRate(outE / GeV / 1e6, TPhotonSpectrum(_fpUniverse->Infrared()->Spectrum(_fPosition.x()/Mpc,_fPosition.y()/Mpc,_fPosition.z()/Mpc)));
  }
  outFile.close();
  std::cout << "Table created " << std::endl;
      throw TCrpErr("STOP Table created");
}

// LM
double TNucleus::SibyllInteraction() {

  TInteractionData *lpInt = _fpUniverse->InteractionData() ;
  TGas* lpGas = _fpUniverse->Gas();

  if (lpGas && (_fpUniverse->Type() != UNIVERSE_ENV1D || _fPosition.x()>0)) {

    if (lpGas->Type() == CLUSTER)
      return 1.0/(lpInt->XSecPPProton(_fEnergy)*lpGas->Density(_fPosition.mag()/Mpc)/centimeter3);
    if (lpGas->Type() == GRIDGAS)
      return 1.0/(lpInt->XSecPPProton(_fEnergy)*lpGas->Density(_fPosition/Mpc)/centimeter3);
  }
  return 1.e99;
}

// LM

static const double convfactor = microbarn/centimeter3*pow(GeV/eV,2);
static const double pm         = 0.93827;
static const double smin       = 1.1646;

double TNucleus::CalcIRLossRate(double lRefEnergy, TPhotonSpectrum infra) {

  double Pp = sqrt(lRefEnergy*lRefEnergy-pm*pm);
  double epsm1 = max(infra.Epsmin(), GeV/eV*0.5*(smin-pm*pm)/(lRefEnergy+Pp));

  //cout << " Epsmax " << infra.Epsmax() << " epsm1 " << epsm1 << endl;

  if (infra.Epsmax() <= epsm1) return 0.0;

  int NE = infra.NE();
  double* infraredfield = new double[2*NE+3]();
  infraredfield[0] = double(NE);
  
  for (int i = 0; i< NE; i++) infraredfield[i+1] = infra.GetE()[i];
  for (int i = 0; i< NE; i++) infraredfield[i+NE+1] = infra.Spectrum()[i];
  
  infraredfield[2*NE+1] = lRefEnergy;
  infraredfield[2*NE+2] = (_fCharge != 0) ? 13 : 14;
  TF1 fti("fti", this, &TNucleus::functs_int_ir, -1e20, 1e20, 2*NE+3, "TNucleus");
  for (int i = 0; i < 2*NE+3; i++) fti.SetParameter(i,infraredfield[i]);
    double result = fti.Integral(log(epsm1), log(infra.Epsmax()), infraredfield, 1e-3);
  
  //cout << " Epsmax " << infra.Epsmax()/EeV << " epsm1 " << epsm1/EeV << endl;
  //  double result = GSLIntegrator(functs_int_ir, log(epsm1), log(infra->Epsmax()), infraredfield);
  
  delete [] infraredfield;
  infraredfield= NULL;

  //  cout << lRefEnergy << " " << infra.Epsmin() << " " << infra.Epsmax() << " " << epsm1 << " " << (result/(8.0*lRefEnergy*Pp)*convfactor)*Mpc << endl;
  return (result/(8.0*lRefEnergy*Pp)*convfactor);
}


double TNucleus::functs_int_ir(double* eps_ln, double* param) {

  //  double* param = (double*)param1;
  TInteractionData* lpInt = _fpUniverse->InteractionData();
  int N = int(param[0]);
  //  std::cout << "functs_int_ir : " << N << std::endl;
  double E0 = param[2*N+1];
  double partid = param[2*N+2];

  double Pp = sqrt(E0*E0-pm*pm);
  double eps=exp(*eps_ln);

  double smax = pm*pm+2.e0*eps*(E0+Pp)*eV/GeV;
  if (smax <= smin) return 0.0;

  //  double* param_xsec = new double[2];
  //param_xsec[0] = double(partid);
  //param_xsec[1] = E0;

  int i = 1;
  double photd = param[N+i];
  double t = 0.0;
  if (param[i] != eps) {
    while (param[i] < eps && i < N+1) i++;
    if (i < N+1) i--;
    t = (eps-param[i])/(param[i+1]-param[i]);
    photd = param[N+i]*(1.0-t) + param[N+i+1]*t;
  }

  int index_E = min(int(log10(E0/lpInt->E_IR(0))/lpInt->dEtabXsecIR()) , 250) ;
  int index_S = min(int((smax-lpInt->S_IR(0))/lpInt->dStabXsecIR()) , 999) ;
  //std::cout << "After dE & dS  index_E : " << index_E << " index_S : " << index_S << std::endl;
  if(index_E < 0)return 0.;
  double result = (partid==13) ? lpInt->IRXsecInt_Proton(index_E, index_S) : lpInt->IRXsecInt_Neutron(index_E, index_S);//GSLIntegrator(func_sigma, smin, smax, param_xsec);

  //  delete [] param_xsec;
  //param_xsec = NULL;
  //std::cout << "end of functs" << std::endl;
  return result*photd/pow(eps,3);
}

// LM
void TNucleus::SibyllPProd() {
  TInteractionData *lpInt = _fpUniverse->InteractionData() ;

  int lNature = 0;
  if (abs(_fCharge) == 0) lNature = 1 ;

  double lEinGeV = _fEnergy* invGeV ;
  double lMomentaList[5][15000] ; 
  int lParticleList[15000] ;
  int lNbParticle, lPartType ;
  int lCheckFlag= 0 ; // to check only one p or n is out!
  cout << lNature << " " << lEinGeV << endl;
  sibyllevent(lNature, lEinGeV, lMomentaList, lParticleList, lNbParticle);

  // important : the incoming momentum is along z for sophia.
  // we do only use the outcoming energy, ie. lMomentaList[3][]
  for (int i=0; i<lNbParticle; i++) {
    lPartType = lParticleList[i] ;
    
    if (lPartType == 13 || lPartType == 14) {

      if (lCheckFlag == 0) {
	lCheckFlag += 1 ;
	if (lPartType == 13) _fCharge = 1 ;
	if (lPartType == 14) _fCharge = 0 ;
	_fEnergy = lMomentaList[3][i] * GeV ;
	_fMomentum.setMag(_fEnergy * inv_c_light) ; 
      } else { // multinucleon production --> need to add a nucleon in the queue
	double lEnergy = lMomentaList[3][i] * GeV ;
	if (lEnergy < _fpUniverse->Sources()->Emin()) continue;
	int lCharge = 0 ;
	if (lPartType == 13) lCharge = 1 ;
	if (lPartType == 14) lCharge = 0 ;
	TVector3D lMomentum = _fMomentum ;
	lMomentum.setMag( lMomentaList[3][i] * GeV * inv_c_light ) ;
	fParticleQueue.push_back(new TNucleus(_fpUniverse, _fPosition, lMomentum, _fInitPosition, 
					      _fInitMomentum, _fTime, 1, lCharge, _fpList1DPhotons, _fInitType
						    )) ;
      }
    } else if (lPartType == 1 || lPartType == 2 || lPartType == 3) {
      // New electromagnetic cascade.
      double lCharge = 0;
      if (lPartType == 2) lCharge = eplus ;
      if (lPartType == 3) lCharge = electron_charge ;
      TVector3D lGammaMomentum = _fMomentum ;
      lGammaMomentum.setMag(lMomentaList[3][i] * GeV * inv_c_light) ;
      if (lpInt->SecPhotonFlag()) 
	if (_fpUniverse->Type() == UNIVERSE_ENV1D) {
	  PARTICLE lPart ;
	  if (lCharge == 0) {
	    lPart = PHOTON ;
	  }else if (lCharge == eplus) {
	    lPart = POSITRON ;
	  } else if (lCharge == electron_charge) {
	    lPart = ELECTRON ;
	  } else throw TCrpErr("Charge error!") ;
	  _fpList1DPhotons->AddPhotons(_fPosition.x(), lGammaMomentum.mag()*c_light, lPart) ;
	} else {
	  string lPart ;
	  if (lCharge == 0) {
	    lPart = "PP_gamma";
	  }else if (lCharge == eplus) {
	    lPart = "PP_e+";
	  } else if (lCharge == electron_charge) {
	    lPart = "PP_e-" ;
	  } else throw TCrpErr("Charge error!") ;
	  fParticleQueue.push_back(new TPhoton(_fPosition, lGammaMomentum, lCharge, _fInitPosition, 
					       _fInitMomentum.getR()*c_light, _fpUniverse, _fInitType, lPart, _fTime));
	}
    } else if (lPartType >= 15 && lPartType <= 18) {
      // New neutrino!
      int lFlavor = 12 ; 
      if (lPartType == 15) lFlavor = 12 ; // nu_e
      if (lPartType == 16) lFlavor = -12 ; // nu_e~
      if (lPartType == 17) lFlavor = 14 ; // nu_mu
      if (lPartType == 18) lFlavor = -14 ; // nu_mu~
      TVector3D lNuMomentum = _fMomentum ;
      lNuMomentum.setMag(lMomentaList[3][i] * GeV * inv_c_light) ;
      if (lpInt->SecNuFlag()) 
	fParticleQueue.push_back(new TNeutrino(lFlavor,_fPosition, lNuMomentum, _fInitPosition, _fInitMomentum.mag()*c_light, _fpUniverse, _fInitType, _fTime)) ;
      
    } else if (lPartType == -13 || lPartType == -14) {
      // Antiproton/neutron --> forget it! (very small flux)
    } else throw TCrpErr("Unknown particle generated by Sibyll.") ;

  }
  if (lCheckFlag == 0) throw TCrpErr("No nucleon generated by Sibyll.") ;

}






double TNucleus::RandPowerLaw(double index, double min, double max, double toto) {
	//check for index -1!
	if ((std::abs(index + 1.0)) < std::numeric_limits<double>::epsilon()) {
		double part1 = log(max);
		double part2 = log(min);
		return exp((part1 - part2) * toto + part2);
	} else {
		double part1 = pow(max, index + 1);
		double part2 = pow(min, index + 1);
		double ex = 1 / (index + 1);
		return pow((part1 - part2) * toto + part2, ex);
	}
}


