/**
   @file    photon.cc
   @author  Eric Armengaud, armengau@in2p3.fr
   @brief   Implementation of the TPhoton class. See the .h file
*/

#include "photon.h"
#include <sstream>

TPhoton::TPhoton(double aPosition,TUniverse* aUniv, char* type) {
  // Trivial constructor for the 1D case

  SetType(PARTICLE_PHOTON);
  _fCharge = 10 ; // special flag...
  _fPosition.set(aPosition,0,0) ;
  _fInitPosition = _fPosition ;
  _fSourcePosition = _fPosition ;
  _fMomentum = TVector3D() ;
  _fEnergy = 0 ;
  _fPairProdFlag = 0 ;
  _fSourceEnergy = 0 ;
  New_dCVector(&_fEnergyGrid, NUM_MAIN_BINS) ;
  New_dCVector(&_fEnergyWidth, NUM_MAIN_BINS);
  SetEnergyBins(MIN_ENERGY_EXP, &_fEnergyGrid, &_fEnergyWidth);
  NewSpectrum(&_fFullSpectrum, NUM_MAIN_BINS);
  _fpUniverse = aUniv;
  _fInitType = 22;
  _fOriginStr=type;
  _fInitTime=0;
    da = _fpUniverse->DistanceArray();
}


TPhoton::TPhoton(TUniverse* aUniv) {
  SetType(PARTICLE_PHOTON);
  _fpUniverse = aUniv;
    da = _fpUniverse->DistanceArray();
  throw TCrpErr("Gamma sources not implemented!") ;
}

TPhoton::TPhoton(TVector3D aPosition, TVector3D aMomentum, double aCharge, TVector3D aSourcePosition, double aSourceEnergy, TUniverse* aUniv, int aInitType, string aOriginStr, double aInitTime) {
  SetType(PARTICLE_PHOTON);
  // Constructor from pion production products
#ifdef DEBUG_OUTPUT
  cout << "TPhoton constructor from pion production" << endl;
#endif
  _fCharge = aCharge ;
  if (aCharge == 0) _fMass = 0 ;
  else if (aCharge == eplus || aCharge == electron_charge) _fMass = electron_mass_c2/(MeV * c_squared) ;
  else throw TCrpErr("charge not ok in photon constructor.") ;
  _fPosition = aPosition ;
  _fInitPosition = aPosition ;
  _fSourcePosition = aSourcePosition ;
  _fSourceEnergy = aSourceEnergy ;
  _fMomentum = aMomentum ;
  _fEnergy = _fMomentum.mag() * c_light ;
  _fPairProdFlag = 0 ;
  _fpUniverse = aUniv;
  _fInitType=aInitType;
  _fOriginStr=aOriginStr;
  _fInitTime=aInitTime;
  // Initialize the spectrum and energy grids for dint
  New_dCVector(&_fEnergyGrid, NUM_MAIN_BINS) ;
  New_dCVector(&_fEnergyWidth, NUM_MAIN_BINS);
  SetEnergyBins(MIN_ENERGY_EXP, &_fEnergyGrid, &_fEnergyWidth);
  NewSpectrum(&_fFullSpectrum, NUM_MAIN_BINS);
    da = _fpUniverse->DistanceArray();
  PARTICLE lPart ;
  if (_fCharge == 0) {
    lPart = PHOTON ;
  } else if (_fCharge == eplus) {
    lPart = POSITRON ;
  } else if (_fCharge == electron_charge) {
    lPart = ELECTRON ;
  } else throw TCrpErr("charge not ok in photon constructor") ;
  this->AddToSpectrum(_fEnergy, lPart) ;

}

TPhoton::TPhoton(double aInjE_hadron, double aDeltaE_hadron, TVector3D aPosition, TVector3D aMomentum, TVector3D aSourcePosition, double aSourceEnergy, TUniverse* aUniv, int aMassNumber, int aInitType, string aOriginStr, double aInitTime) {
  SetType(PARTICLE_PHOTON);
  // Constructor from proton pair production

  _fCharge = 0 ; // e+e- pair = 0 charge (convention)
  _fPosition = aPosition ;
  _fInitPosition = aPosition ;
  _fDeltaE_hadron = aDeltaE_hadron ;
  _fInjE_hadron = aInjE_hadron ;
  _fEnergy = aDeltaE_hadron ;
  _fMomentum = aMomentum ;
  _fSourcePosition = aSourcePosition ;
  _fSourceEnergy = aSourceEnergy ;
  _fPairProdFlag = 1 ;
  _fpUniverse = aUniv;
  _fInitType=aInitType;
  _fOriginStr=aOriginStr;
  _fInitTime=aInitTime;
    da = _fpUniverse->DistanceArray();
  // Initialize the spectrum and energy grids for dint
  New_dCVector(&_fEnergyGrid, NUM_MAIN_BINS) ;
  New_dCVector(&_fEnergyWidth, NUM_MAIN_BINS);
  SetEnergyBins(MIN_ENERGY_EXP, &_fEnergyGrid, &_fEnergyWidth);
  NewSpectrum(&_fFullSpectrum, NUM_MAIN_BINS);
  this->AddToSpectrum(_fInjE_hadron, _fDeltaE_hadron, aMassNumber) ;

}

TPhoton::TPhoton(TUniverse* aUniv, double aInjE_hadron, double aDeltaE_hadron, double aTimeStep, TVector3D aPosition, TVector3D aMomentum, TVector3D aSourcePosition, double aSourceEnergy, double aRedshift, int aMassNumber, int aInitType, string aOriginStr, double aInitTime) {
  SetType(PARTICLE_PHOTON);
  // Constructor from proton pair production (Kelner parametrization)

  _fCharge = 0 ; // e+e- pair = 0 charge (convention)
  _fPosition = aPosition ;
  _fInitPosition = aPosition ;
  _fDeltaE_hadron = aDeltaE_hadron ;
  _fInjE_hadron = aInjE_hadron ;
  _fEnergy = aDeltaE_hadron ;
  _fMomentum = aMomentum ;
  _fSourcePosition = aSourcePosition ;
  _fSourceEnergy = aSourceEnergy ;
  _fPairProdFlag = 1 ;
  _fpUniverse = aUniv;
  _fInitType=aInitType;
  _fOriginStr=aOriginStr;
  _fInitTime=aInitTime;
    da = _fpUniverse->DistanceArray();
  // Initialize the spectrum and energy grids for dint
  New_dCVector(&_fEnergyGrid, NUM_MAIN_BINS) ;
  New_dCVector(&_fEnergyWidth, NUM_MAIN_BINS);
  SetEnergyBins(MIN_ENERGY_EXP, &_fEnergyGrid, &_fEnergyWidth);
  NewSpectrum(&_fFullSpectrum, NUM_MAIN_BINS);
  this->AddPairProdSpectrum(_fInjE_hadron, aTimeStep, aRedshift, aUniv) ;

}


void TPhoton::AddToSpectrum(double aInjEnergy, PARTICLE aPart) {
  // Add a single particle (from pion production or direct photon injection)
  // --> monochromatic injection

  int num_main_bins = _fEnergyGrid.dimension ;
  if ((_fEnergyWidth.dimension != num_main_bins) ||
      (_fFullSpectrum.numberOfMainBins != num_main_bins))
    throw TCrpErr("TPhoton spectra : inconsistent dimensions");
  
  double criticalEnergy = aInjEnergy/(eV * ELECTRON_MASS) ; // units of dint
  int maxBin = (int)((log10(criticalEnergy*ELECTRON_MASS) -
		    MAX_ENERGY_EXP)*BINS_PER_DECADE + num_main_bins);
  (_fFullSpectrum.spectrum)[aPart][maxBin] += 1.; // We add (and not set) 1 particle in the bin
    
}

void TPhoton::AddToSpectrum(double aEmin, double aEmax, double aEint, double aAlpha) {
  // Add to the spectrum a power-law distribution of photons, of
  // spectral index aAlpha, between aEmin and aEmax, with total energy aEint

  if (aEint == 0.e0) throw TCrpErr("aEint = 0 in TPhoton");
  double sum=0.;
  int num_main_bins = _fEnergyGrid.dimension ;
  if ((_fEnergyWidth.dimension != num_main_bins) ||
      (_fFullSpectrum.numberOfMainBins != num_main_bins))
    throw TCrpErr("TPhoton spectra : inconsistent dimensions");
  double lEint = aEint/(eV*ELECTRON_MASS) ; // conversion to dint system
  double lEmin = aEmin/(eV*ELECTRON_MASS) ;
  double lEmax = aEmax/(eV*ELECTRON_MASS) ;
  Spectrum lSpectrum ; // local spectrum to build and then add to _fFullSpectrum
  NewSpectrum(&lSpectrum, NUM_MAIN_BINS); // la il faut aussi NUM_SPECIES.

  for (int i=0; i<num_main_bins; i++) {
    if (_fEnergyGrid.vector[i] < lEmax && _fEnergyGrid.vector[i] > lEmin) {
      (lSpectrum.spectrum)[PHOTON][i] = pow(_fEnergyGrid.vector[i],-aAlpha)*
	  (_fEnergyWidth.vector)[i];
	sum += (lSpectrum.spectrum)[PHOTON][i]*(_fEnergyGrid.vector)[i];
    }
  }
  sum = lEint/sum ;
  for (int i=0; i<num_main_bins; i++) {
    (lSpectrum.spectrum)[PHOTON][i] *= sum;
    (_fFullSpectrum.spectrum)[PHOTON][i] += (lSpectrum.spectrum)[PHOTON][i] ;
  }
  DeleteSpectrum(&lSpectrum);
}

void TPhoton::AddToSpectrum(double aHInjEnergy, double aDeltaE_hadron, int aMassNumber) {
  // Add to the spectrum a pair production spectrum
  // with a power law of index PAIRPROD_SPECTRALINDEX
  
  if (aDeltaE_hadron == 0.e0) throw TCrpErr("aDeltaE_Hadron = 0 in TPhoton");
  double sum=0.;
  int num_main_bins = _fEnergyGrid.dimension ;
  if ((_fEnergyWidth.dimension != num_main_bins) ||
      (_fFullSpectrum.numberOfMainBins != num_main_bins))
    throw TCrpErr("TPhoton spectra : inconsistent dimensions");

  double criticalEnergy = aHInjEnergy/(eV*ELECTRON_MASS) ; // conversion to dint system

  criticalEnergy = 4.*(aHInjEnergy/eV)*(aHInjEnergy/eV)*1.e-3/
    (4.*(aHInjEnergy/eV)*1.e-3+(8.8e17 * aMassNumber * aMassNumber))/ELECTRON_MASS; // threshold energy (Gunter)
  // 1.e-3 = epsilon/eV ; 9e17 = (mp/eV)^2
 double critEnergy_inf = (3.3e13)/ELECTRON_MASS; // low-energy threshold (ADDED), dint unit
  double lDeltaE_hadron = aDeltaE_hadron/(eV*ELECTRON_MASS) ; // dint units
  Spectrum lSpectrum ; // local spectrum to build and then add to _fFullSpectrum
  NewSpectrum(&lSpectrum, NUM_MAIN_BINS); // la il faut aussi NUM_SPECIES.

  for (int i=0; i<num_main_bins; i++) {
    //    if (_fEnergyGrid.vector[i] < criticalEnergy) {
   if (_fEnergyGrid.vector[i] < criticalEnergy && _fEnergyGrid.vector[i] > critEnergy_inf) {
      (lSpectrum.spectrum)[ELECTRON][i] = pow(_fEnergyGrid.vector[i],PAIRPROD_SPECTRALINDEX)*
	  (_fEnergyWidth.vector)[i];
	sum += (lSpectrum.spectrum)[ELECTRON][i]*(_fEnergyGrid.vector)[i];
    }
  }
  if((sum != 0 )) {
    
    sum = lDeltaE_hadron/sum/2.; // factor 2 = repartition of energy in e+ and e-
    for (int i=0; i<num_main_bins; i++) {
      (lSpectrum.spectrum)[ELECTRON][i] *= sum;
      (lSpectrum.spectrum)[POSITRON][i] = (lSpectrum.spectrum)[ELECTRON][i];
      (_fFullSpectrum.spectrum)[ELECTRON][i] += (lSpectrum.spectrum)[ELECTRON][i] ;
      (_fFullSpectrum.spectrum)[POSITRON][i] += (lSpectrum.spectrum)[POSITRON][i] ;
    }
  }
  DeleteSpectrum(&lSpectrum);
  
  //TESTING
  // for(int i=0; i<num_main_bins; i++){
  //   if (!(
  // 	  (_fFullSpectrum.spectrum)[ELECTRON][i]==(_fFullSpectrum.spectrum)[ELECTRON][i] ||
  // 	  (_fFullSpectrum.spectrum)[POSITRON][i]==(_fFullSpectrum.spectrum)[POSITRON][i] )){ 
      
  //     std::stringstream out;
  //     out <<  i << " aHInjEnergy " << aHInjEnergy << " aDeltaE_hadron " << aDeltaE_hadron << " aMassNumber " << aMassNumber ;
  //     throw TCrpErr("AddToSpectrum nan i " + out.str());
  //   }
  // }
}

void TPhoton::AddPairProdSpectrum(double aEp, double aTimeStep, double aRedshift, const TUniverse* aUniv) {

  TInteractionData *lpInt = aUniv->InteractionData() ;
  double lEp_eV=aEp/eV;
  double lTimeStep_sec=aTimeStep/second ; // units for the spectrum
  double lPairProdSpec=0, lEe=0;
  lEp_eV *= (1+aRedshift) ; // z effect

  for (int i=0; i<NUM_MAIN_BINS; i++) {

    lEe = _fEnergyGrid.vector[i];
    int i_redshift=i; // redshift formula
    // dN/dE(Ep,Ee,z) = (1+z)^4 * dN/dE(Ep*(1+z),Ee*(1+z),0)
    if (aRedshift != 0) {
      lEe *= (1+aRedshift);
      while (lEe > (_fEnergyGrid.vector)[i_redshift]+(_fEnergyWidth.vector)[i_redshift]/2.) {
	i_redshift++;
	if (i_redshift == NUM_MAIN_BINS-1) break;
      }
    }

    // Use the Kelner formula only in its range of validity.
    if (lEp_eV < PAIRPRODCUT_EE_EP*lEe*ELECTRON_MASS) lPairProdSpec=0;
    else {
      lPairProdSpec = pow(1+aRedshift,4)*lpInt->PairProdSpec(lEp_eV,i_redshift) ; // eV/sec
      lPairProdSpec *= lTimeStep_sec;
      lPairProdSpec *= ((_fEnergyWidth.vector)[i]*ELECTRON_MASS) ; // valeurs en eV
      lPairProdSpec /= pow(_fEnergyGrid.vector[i]*ELECTRON_MASS,2);
    }
    (_fFullSpectrum.spectrum)[ELECTRON][i] += lPairProdSpec;
    (_fFullSpectrum.spectrum)[POSITRON][i] += lPairProdSpec;
  }

}


QUEUE<TParticle*>* TPhoton::Propagate(TBasicParam* aBasic) {
  //#define GAM_TEST
#ifdef GAM_TEST
  _fInitPosition.set(3 * Mpc, 3 * Mpc, 3 * Mpc);
  _fPosition=_fInitPosition;
  _fMomentum.setRThetaPhi( _fEnergy* inv_c_light, 3.141592653 * 0 * ( 0.1475846), 3.141592653 * 0/2. ) ;
#endif

#ifdef DEBUG_OUTPUT
  std::cout << "Start of photon propagation" << std::endl; 
#endif
vector<double> lDistanceArray = this->CheckDetection(aBasic) ;
 TVector3D _fPosOld=_fPosition;
 TVector3D _fInitPosOld=_fInitPosition;

  for(int i=0; i< lDistanceArray.size() ; i++){
    _fShowerPropDistance=lDistanceArray[i]; //Reuse the old code: just one detection per photon
    _fPosition=_fPosOld;
    _fInitPosition=_fInitPosOld;
    _fShowerSpectrum.clear();
    //std::cout << "NO Shower development!!!" << std::endl; 
    this->DevelopShower() ;
    //_fPosition += _fMomentum.unit() * _fShowerPropDistance ; // important correction to observer
    this->WriteShower(aBasic) ;
  }
#ifdef GAM_TEST
  exit(-1);
#endif

#ifdef DEBUG_OUTPUT
  std::cout << "End of photon propagation "<< std::endl;
#endif
  return &fParticleQueue ;
}

vector<double> TPhoton::CheckDetection(TBasicParam* aBasic) {

vector<double> lReturn;
  lReturn.clear();
  if ( _fpUniverse->Observers()->Type() == OBSERVER_POINT ) {
    _fShowerPropDistance = _fPosition.x() ;
    lReturn.push_back(_fPosition.x());

  } else if ( _fpUniverse->Observers()->Type() == OBSERVER_SMALLSPHERE ) {

    TObservers* lObs = _fpUniverse->Observers() ;
    double lDistObs;
    vector<TVector3D> lTranslators;
    TVector3D lDum ;
    TVector3D lPos=_fPosition; 
    double dt, lt; 
    TVector3D lBox; 
    TVector3D beta=_fMomentum / _fEnergy * c_light;
    TVector3D lIndex;
    lt=0;

    //Check where 
    lBox.setX( (_fMomentum / _fEnergy).x() < 0? _fpUniverse->Xmin() : _fpUniverse->Xmax() );
    lBox.setY( (_fMomentum / _fEnergy).y() < 0? _fpUniverse->Ymin() : _fpUniverse->Ymax() );
    lBox.setZ( (_fMomentum / _fEnergy).z() < 0? _fpUniverse->Zmin() : _fpUniverse->Zmax() );

    lDum.set(0,0,0);
    lTranslators.push_back( lDum);
    lIndex.set(0, 0, 0); 
    int counter=0;

//     if (PHOTON_DETECTOR_CROWN) {
//       for (int ax=-1; ax<2; ax++) {
// 	for (int ay=-1; ay<2; ay++) {
// 	  for (int az=-1; az<2; az++) {
// 	    lDum.set(ax*(_fpUniverse->Xmax()-_fpUniverse->Xmin()),
// 		     ay*(_fpUniverse->Ymax()-_fpUniverse->Ymin()),
// 		     az*(_fpUniverse->Zmax()-_fpUniverse->Zmin())) ;
// 	    lTranslators.push_back(lDum) ;
// 	  }
// 	}
//       }
//     } else {
//       lDum.set(0,0,0) ;
//       lTranslators.push_back(lDum) ;
//     }

    while(lt < aBasic->Tmax() - _fInitTime){
#ifdef GAM_TEST
      if(++counter > 10) break;
      std::cout << "lt " << lt / Mpc * c_light << "  aBasic->Tmax() " << aBasic->Tmax() / Mpc * c_light << " _fInitTime " << _fInitTime / Mpc * c_light << std::endl;
      std::cout << " aBasic->Tmax() - _fInitTime " <<  (aBasic->Tmax() - _fInitTime)/Mpc * c_light << std::endl;
#endif
      dt=           (lBox.x()-lPos.x())/beta.x() / c_light;
      dt= min(dt, ( (lBox.y()-lPos.y())/beta.y() /c_light));
      dt= min(dt, ( (lBox.z()-lPos.z())/beta.z() /c_light));
#ifdef GAM_TEST
      std::cout << "beta " << beta << " dt " << dt / Mpc * c_light << "  (lBox.z()-lPos.z())/beta.z() " <<  (lBox.z()-lPos.z())/beta.z() / Mpc * c_light  << std::endl;
      std::cout << " lBox.z() " << lBox.z()/Mpc << " lPos.z() " << lPos.z()/Mpc << std::endl;
      std::cout << " lPos " << lPos / Mpc << std::endl;
#endif

      if(dt < 1e-4 *Mpc / c_light) dt=1e-4 * Mpc / c_light;
      lPos+=dt * beta * c_light * (1 + 1e-5);
#ifdef GAM_TEST
      std::cout << " lPos " << lPos / Mpc << std::endl;
#endif
      lt+=dt;

      int lflag=0; 
      if(lPos.x() > _fpUniverse->Xmax()) {
	lPos.setX(lPos.x() -  (_fpUniverse->Xmax()-_fpUniverse->Xmin())) ;
	lIndex.setX(lIndex.x()+1);
	lflag++	;
      } else if(lPos.x() < _fpUniverse->Xmin()) {
	lPos.setX(lPos.x() +  (_fpUniverse->Xmax()-_fpUniverse->Xmin())) ;
	lIndex.setX(lIndex.x()- 1);
	lflag++	;
      }

      if(lPos.y() > _fpUniverse->Ymax()) {
	lPos.setY(lPos.y() -  (_fpUniverse->Ymax()-_fpUniverse->Ymin())) ;
	lIndex.setY(lIndex.y()+1);
	lflag++	;
      } else if(lPos.y() < _fpUniverse->Ymin()) {
	lPos.setY(lPos.y() +  (_fpUniverse->Ymax()-_fpUniverse->Ymin())) ;
	lIndex.setY(lIndex.y()- 1);
	lflag++;	
      }

      if(lPos.z() > _fpUniverse->Zmax()) {
	lPos.setZ(lPos.z() -  (_fpUniverse->Zmax()-_fpUniverse->Zmin())) ;
	lIndex.setZ(lIndex.z()+1);
	lflag++;	
      } else if(lPos.z() < _fpUniverse->Zmin()) {
	lPos.setZ(lPos.z() +  (_fpUniverse->Zmax()-_fpUniverse->Zmin())) ;
	lIndex.setZ(lIndex.z()- 1);
	lflag++; 	
      }

      if(lflag == 0) std::cout << " lflag == 0 " << std::endl;
      lDum.setX(lIndex.x() * (_fpUniverse->Xmax()-_fpUniverse->Xmin()));
      lDum.setY(lIndex.y() *(_fpUniverse->Ymax()-_fpUniverse->Ymin()));
      lDum.setZ(lIndex.z() *(_fpUniverse->Zmax()-_fpUniverse->Zmin()));

      long storA, storB;
      try{
	storA=lTranslators.size();
	storB=lTranslators.max_size();
	if(lt < aBasic->Tmax() - _fInitTime) lTranslators.push_back(lDum) ;
      }    catch (exception& e)
	{
	  std::cout << " 4rd blck " <<  e.what() << std::endl;
	  std::cout << " lTranslators.size() " << lTranslators.size() << " lTranslators.max_size() " << lTranslators.max_size() << std::endl;
	  std::cout << " lTranslators.max_size() - lTranslators.size() " << lTranslators.max_size() - lTranslators.size() << std::endl; 
	  std::cout << " - storA=lTranslators.size(); + storB=lTranslators.max_size(); " << storB - storA << std::endl;
	  exit(-1);
	}  

    }


#ifdef GAM_TEST
    std::cout << "lt " << lt /Mpc * c_light << std::endl;
    std::cout << "lTranslators.size() : " << lTranslators.size() << std::endl;
#endif


    
    for (unsigned int k_trans=0; k_trans<lTranslators.size(); k_trans++) {
#ifdef GAM_TEST
      std::cout << lTranslators[k_trans].x()/ Mpc << " " << lTranslators[k_trans].y()  / Mpc << " " << lTranslators[k_trans].z() / Mpc << std::endl;
#endif
      for (int j=0; j<lObs->Nb(); j++) {
	TVector3D lObsPart = _fPosition - (lObs->Positions(j)+lTranslators[k_trans]) ;
	lDistObs = lObsPart.mag() - lObs->Radius() ;
	if (lDistObs > 0 ) {
	  double lXP, lDiscrim ;
	  lXP = lObsPart.dot(_fMomentum.unit()) ;
	  lDiscrim = pow(lXP,2) - lObsPart.mag2() + pow(lObs->Radius(),2) ;
	  if ( lDiscrim >= 0 && lXP <= 0) { // Direction defined by momentum hits sphere: photon hits the detector.
	    //lReturn = 1 ;
	    //	    _fShowerPropDistance = min(_fShowerPropDistance,lDistObs);
	    if((-lXP-sqrt(lDiscrim)) < (aBasic->Tmax() - _fInitTime)*c_light) lReturn.push_back(-lXP-sqrt(lDiscrim));
	    // hence we select the nearest observer for shower detection.
	  }
	}
      }
    }
  } else if ( _fpUniverse->Observers()->Type() == OBSERVER_LARGESPHERE ) {
    TObservers* lObs = _fpUniverse->Observers() ;
    double lDistObs ;
    int j=0 ; // Only one observer
    lDistObs = (_fPosition-_fSourcePosition).mag() - lObs->Radii(j) ;
    if (lDistObs < 0) {//Photon inside of the sphere
      double lCosTheta = ((_fPosition-_fSourcePosition).unit()).dot(_fMomentum.unit()) ;
      double lSP = (_fPosition-_fSourcePosition).mag() ;
      double lDist = sqrt(pow(lObs->Radii(j),2) - lSP*lSP*(1-lCosTheta*lCosTheta)) - lSP*lCosTheta ;
      lReturn.push_back(lDist);
      _fShowerPropDistance = lDist ;
    }

  } else if ( _fpUniverse->Observers()->Type() == OBSERVER_NO ) {
    lReturn.clear();

  } else throw TCrpErr("Error : CheckDetection for photons not implemented in this configuration" );

#ifdef GAM_TEST
  for(int k=0; k< lReturn.size() ; k++) std::cout<< lReturn[k] / Mpc << std::endl;
#endif
  return lReturn ;
}

void TPhoton::DevelopShower() {
  string lDirTable = _fpUniverse->InteractionData()->ShowerTableDir();

  // 1) Build the B_array which is necessary.
  TMagField* lpField = _fpUniverse->MagField() ;
  dCVector lBPerp ;

  if (lpField->Type() == MAGFIELD_NO) {
    New_dCVector(&lBPerp, 5);
    for (int i=0; i<5; i++) (lBPerp.vector)[i] = 0.;

  } else if (lpField->Type() == MAGFIELD_1D) {
    int lN = (int)(_fShowerPropDistance/lpField->Stepsize());
    if (lN == 0) lN = 1 ;
    New_dCVector(&lBPerp,lN) ;
    TVector3D lPosition ;
    for (int i=0; i<lN; i++) {
      lPosition=_fPosition+_fMomentum.unit()*_fShowerPropDistance*(i+0.5)/lN;
      (lBPerp.vector)[i]=(lpField->getField(lPosition)).perp(_fMomentum)/gauss;
    }

  } else if (lpField->Type() == MAGFIELD_UNIFORM) {
    New_dCVector(&lBPerp, 5);
    double lBPerpGauss = 
      (lpField->getField(_fPosition)).perp(_fMomentum)/gauss;
    for (int i=0; i<5; i++) (lBPerp.vector)[i] = lBPerpGauss ;

  } else if (lpField->Type() == MAGFIELD_LSS || 
	     lpField->Type() == MAGFIELD_KOLMOGOROFF) {
    int lN = (int)(_fShowerPropDistance/lpField->Stepsize());
    New_dCVector(&lBPerp,lN) ;
    TVector3D lPosition ;
    for (int i=0; i<lN; i++) {
      lPosition=_fPosition+_fMomentum.unit()*_fShowerPropDistance*(i+0.5)/lN;
      (lBPerp.vector)[i]=(lpField->getField(lPosition)).perp(_fMomentum)/gauss;
    }
  } else throw TCrpErr("mag field not ok with shower of photons...") ;

  // 2) Call prop_second routine
  int lIRFlag = 0 ;
  int lRadioFlag = 0 ; // Flags corresponding to the functions of prop_second
  if (_fpUniverse->Infrared()->Type() == IR_UNIFORM_HIGH) lIRFlag = 0 ;
  if (_fpUniverse->Infrared()->Type() == IR_UNIFORM_LOW) lIRFlag = 1 ;
  if (_fpUniverse->Infrared()->Type() == IR_UNIFORM_PRIMACK) lIRFlag = 2 ;
  //  if (lIRFlag != 2) cerr << "Warning : Low/high IR for dint not available otherwise.." << endl;
  if (_fpUniverse->RadioBackground() == "High") lRadioFlag = 0 ;
  if (_fpUniverse->RadioBackground() == "Med") lRadioFlag = 1 ;
  if (_fpUniverse->RadioBackground() == "Obs") lRadioFlag = 2 ;
  if (_fpUniverse->RadioBackground() == "Null") lRadioFlag = 3 ;
  // Cosmological parameters
  double lH0 = DEFAULT_H_0_KM_S_MPC ;
  double lOmegaLambda = DEFAULT_OMEGA_LAMBDA ;
  double lOmegaM = DEFAULT_OMEGA_M ;
  double lZmax = _fpUniverse->Infrared()->Zmax() ;
  if (_fpUniverse->Type() == UNIVERSE_ENV1D) {
    lH0 = _fpUniverse->H0() * second * Mpc / kilometer ;
    lOmegaLambda = _fpUniverse->OmegaLambda() ;
    lOmegaM = _fpUniverse->OmegaM() ;
  }
  // Flag to cut the e+/- cascade due to deflections
  double lCutcascade_Magfield = _fpUniverse->InteractionData()->CutcascadeFlag() ;

  Spectrum lSpectrumOut ; // full spectrum after propagation
  NewSpectrum(&lSpectrumOut, NUM_MAIN_BINS) ;
  prop_second(_fShowerPropDistance* invMpc,
	      &lBPerp,&_fEnergyGrid,
	      &_fEnergyWidth,&_fFullSpectrum, 
	      &lSpectrumOut, lDirTable,
	      lIRFlag, lZmax, lRadioFlag,
	      lH0, lOmegaM, lOmegaLambda,
	      lCutcascade_Magfield) ;

  // 3) Build arrays of spectrum: photon = particle of type 0
  for (int j=0; j<lSpectrumOut.numberOfMainBins; j++) {
    _fShowerSpectrum.push_back((lSpectrumOut.spectrum)[0][j]) ; // spectrum: mean number of particles per energy bin
    //_fShowerSpectrum.push_back((lSpectrumOut.spectrum)[1][j]) ; // spectrum: mean number of particles per energy bin
    //std::cout << "WARNING: DINT photons are electron!!!" << std::endl;
  }

  DeleteSpectrum(&lSpectrumOut) ;
  Delete_dCVector(&lBPerp) ;
  _fPosition += _fMomentum.unit() * _fShowerPropDistance ; // important correction to observer

}

void TPhoton::WriteShower(TBasicParam* aBasic) const {

  TOutputData* lOutput = aBasic->OutputData() ;
  string lOrigin=(_fOriginStr + "          ").substr(0, 10);//lOrign should be exactly ten chars.  
  lOrigin[10]='\0';
  if (_fpUniverse->Type() == UNIVERSE_ENV1D) {
    lOutput->Add1DShower(lOrigin,_fSourcePosition.mag(),_fInitPosition.mag(),_fEnergy,
			 _fSourceEnergy,_fShowerSpectrum) ;
  } else if (_fpUniverse->Type() == UNIVERSE_LARGESCALE) {
    lOutput->Add3DShower(lOrigin,_fSourcePosition,_fInitPosition,_fPosition,
			 _fMomentum,_fSourceEnergy,_fShowerSpectrum, _fInitType) ;
  } else throw TCrpErr( "Detection not implemented in this environment" );
  
}

TPhoton::~TPhoton() {
  _fShowerEnergy.clear() ;
  _fShowerSpectrum.clear() ;
  Delete_dCVector(&_fEnergyGrid) ;
  Delete_dCVector(&_fEnergyWidth) ;
  DeleteSpectrum(&_fFullSpectrum) ;

}
