
/**
   @file    list1dphotons.cc
   @author  Eric Armengaud, armengau@in2p3.fr
   @brief   Implementation of the TList1DPhotons class. See the .h file
*/

#include "list1dphotons.h"

TList1DPhotons::TList1DPhotons(TUniverse* aUniv, char* type) {

  if (aUniv->Type() == UNIVERSE_ENV1D && (aUniv->InteractionData()->SecPairProdPhotonFlag()
					  || aUniv->InteractionData()->SecPhotonFlag()
					  || aUniv->Sources()->IsPhotonSources())) {

    _fStepShowers = DEFAULT_STEP_1DPHOTONSHOWERS_MPC*Mpc ;
      
    if (aUniv->InteractionData()->SecPairProdPhotonFlag())
      _fStepShowers = aUniv->InteractionData()->InteractionTimeStep()*c_light ;

    bool lChooseGrid =  USE_GRID_POSITIONS_1DPHOTONS ;
    double lListPos[] = {1.00000, 3.00000, 5.00000, 7.00000, 9.00000, 11.0000,
			 13.0000, 15.0000, 17.0000, 19.0000, 21.0000, 23.0761,
			 25.1191, 27.3431, 29.7639, 32.3990, 35.2674, 38.3898,
			 41.7887, 45.4884, 49.5157, 53.8995, 58.6715, 63.8660,
			 69.5203, 75.6753, 82.3751, 89.6682, 97.6069, 106.249,
			 115.655, 125.895, 137.041, 149.174, 162.381, 176.757,
			 192.406, 209.440, 227.983, 248.168, 270.139, 294.056,
			 320.090, 348.429, 379.277, 412.856, 449.408, 489.196,
			 532.507, 579.652, 630.971, 686.834, 747.643, 813.835,
			 885.888, 964.320, 1049.69, 1142.63, 1243.79, 1353.91,
			 1473.78, 1604.26, 1746.29, 1900.90, 2069.19, 2252.39,
			 2451.80, 2668.87, 2905.16, 3162.36, 3442.35, 3747.11,
			 4078.86} ;

    if (lChooseGrid) {
      for (int i=0; i<73; i++) {
	if (lListPos[i] <= aUniv->Xmax()/Mpc) {
	   _fListInitPositions.push_back(lListPos[i]*Mpc) ;
	   _fListPhotons.push_back(new TPhoton(_fListInitPositions.at(i), aUniv, type)) ;
	}
      }
      _fNbPhotons = _fListPhotons.size() ;
      _fStepShowers = 0 ;
    } else {
      _fNbPhotons = (int)(aUniv->Xmax()/_fStepShowers)+1 ;
      for (int i=0; i<_fNbPhotons; i++) {
	_fListInitPositions.push_back((i+0.5)*_fStepShowers) ;
	_fListPhotons.push_back(new TPhoton(_fListInitPositions.at(i), aUniv, type)) ;
      }
    }

  } else {
    // Nothing
  }

}

void TList1DPhotons::AddPhoton(const TUniverse* aUniv) {

  // 1 ) Initial position (and energy)
  TVector3D lPosition ;
  double lEnergy ;
  TSources *sources = aUniv->Sources() ;
  if (sources->Type() == SOURCE_DISCRETE) {
    // Pick up a source from the list
    int i = RandFlat::shootInt(sources->Nb()) ;
    lPosition = sources->Positions(i) ;
    
    if ( sources->SpectrumFlag() == 1) { // monochromatic spectrum
      lEnergy = sources->Ecut() ;
    } else if (sources->SpectrumFlag() == 2) { // spectral index
      lEnergy = -1. ;
    }
  } else if (sources->Type() == SOURCE_CONTINUOUS) {
    lPosition = sources->GetPosition() ;
    
    if ( sources->SpectrumFlag() == 1) { // monochr
      lEnergy = sources->Ecut() ;
    } else if (sources->SpectrumFlag() == 2) { // spectral index
      lEnergy = -1. ;
    }
  } else throw TCrpErr("Unknown source type in 1D photon generator." );
  
  // 2) Filling the photon list
  if (sources->SpectrumFlag() == 1) {
    PARTICLE lPart = PHOTON ;
    this->AddPhotons(lPosition.x(), lEnergy, lPart) ;
  } else {
    double lIntensity = 1. ;
    double lEmin = aUniv->Sources()->Emin() ;
    double lEmax = aUniv->Sources()->Ecut() ;
    double lAlpha = aUniv->Sources()->Alpha() ;
    this->AddPhotons(lPosition.x(), lEmin, lEmax, lAlpha, lIntensity) ;
  }

}

void TList1DPhotons::AddPhotons(double aPositionX, double aEmin, double aEmax, double aAlpha, double aIntensity) {
  // Add power-law spectrum

  int lBin = 0 ;
  if (_fStepShowers != 0) {
    lBin = (int)(aPositionX / _fStepShowers) ;
  } else {
    lBin = _fNbPhotons -1 ;
    while ( lBin!=0 && _fListInitPositions.at(lBin) > aPositionX) lBin -= 1 ;
  }
  (_fListPhotons.at(lBin))->AddToSpectrum(aEmin, aEmax, aIntensity, aAlpha) ;

}

void TList1DPhotons::AddPhotons(double aPositionX, double aHInjEnergy, double aDeltaE_hadron, int aMassNumber) {
  // Add spectrum from PPP

  int lBin = 0 ;
  if (_fStepShowers != 0) {
    lBin = (int)(aPositionX / _fStepShowers) ;
  } else {
    lBin = _fNbPhotons -1 ;
    while ( lBin!=0 && _fListInitPositions.at(lBin) > aPositionX) lBin -= 1 ;
  }
  //std::cout << "AddPairProdSpectrum  aHInjEnergy " << aHInjEnergy << " aDeltaE_hadron  " << aDeltaE_hadron << " aMassNumber  " << aMassNumber << std::endl;
  (_fListPhotons.at(lBin))->AddToSpectrum(aHInjEnergy, aDeltaE_hadron, aMassNumber) ;
}

void TList1DPhotons::AddPairProdPhotons(double aPositionX, double aHInjEnergy, double aTimeStep, const TUniverse* aUniv, double aRedshift) {
  // Add spectrum from PPP : using Kelner formula

  int lBin = 0 ;
  if (_fStepShowers != 0) {
    lBin = (int)(aPositionX / _fStepShowers) ;
  } else {
    lBin = _fNbPhotons -1 ;
    while ( lBin!=0 && _fListInitPositions.at(lBin) > aPositionX) lBin -= 1 ;
  } 
  (_fListPhotons.at(lBin))->AddPairProdSpectrum(aHInjEnergy, aTimeStep, aRedshift, aUniv) ;

}

void TList1DPhotons::AddPhotons(double aPositionX, double aInjEnergy, PARTICLE aPart) {
  // Add a single particle

  int lBin = 0 ;
  if (_fStepShowers != 0) {
    lBin = (int)(aPositionX / _fStepShowers) ;
  } else {
    lBin = _fNbPhotons -1 ;
    while ( lBin!=0 && _fListInitPositions.at(lBin) > aPositionX) lBin -= 1 ;
  }
  (_fListPhotons.at(lBin))->AddToSpectrum(aInjEnergy, aPart) ;

}

void TList1DPhotons::Propagate(const TUniverse* aUniv, TBasicParam* aBasic) {

  if (aUniv->Type() == UNIVERSE_ENV1D && (aUniv->InteractionData()->SecPairProdPhotonFlag()
					  || aUniv->InteractionData()->SecPhotonFlag()
					  || aUniv->Sources()->IsPhotonSources())) {

    for (int i=0; i<_fNbPhotons; i++) {
      TPhoton* lpPhoton = _fListPhotons.at(i) ;
      Spectrum lSpectrum = lpPhoton->FullSpectrum() ;
      dCVector lEnergyGrid = lpPhoton->EnergyGrid() ;
      double lEnergy = GetEnergy(&lSpectrum, &lEnergyGrid) ; 
      if (lEnergy != 0) QUEUE<TParticle*>* lDum = lpPhoton->Propagate(aBasic) ;
    }
  } else {
    // Nothing
  }

}

TList1DPhotons::~TList1DPhotons() {
  
  while (_fListPhotons.size()) {
    delete _fListPhotons.back() ;
    _fListPhotons.pop_back() ;
  }
  
}
