

/**
   @file   discretesources.cc
   @author Eric Armengaud, armengau@in2p3.fr
   @brief  Implementation of the TDiscreteSources class. See the .h file
   
*/

//changed by Nils Nierstenhoefer (NN) in 09/2009 to define sources which emit paticles with a fixed mass.

#include "discretesources.h"

TDiscreteSources::TDiscreteSources(const char* aFileName) : TXmlParam(aFileName) {
  SetType(SOURCE_DISCRETE);
  TiXmlElement* lpXmlSources = XmlExtract().GetElement("Sources") ;
  TiXmlElement* lpXmlUniverse = XmlExtract().GetElement("Environment") ;
  string lUnivType = lpXmlUniverse->Attribute("type") ;
  string dum;
  // Read minimal energy
    _fEmin = XmlExtract().GetDouble("MinEnergy_EeV") ;

    if (_fEmin <= 0) throw TXmlErr("Error while parsing minimum energy value.");
    _fEmin *= EeV ;
  
  // Check if this is a photon or nucleus source
  TiXmlElement* lpXmlParticleType = lpXmlSources->FirstChildElement("Particles") ;
  if (lpXmlParticleType) {
    _fInitialMassNumber.push_back(0);
    _fInitialChargeNumber.push_back(0);
    string lParticleType = lpXmlParticleType->Attribute("type") ;
    if (lParticleType == "Photons") {
      _fIsPhotonSources = 1 ;
    } else if (lParticleType == "Protons") {
      //TODO: delete this and use Nucleus instead.
      throw TXmlErr("Particle type \"Protons\" has been choosen. This is obsolete now: Use Nuclei with MassNumber=1 and ChargeNumber=1." );
      _fIsPhotonSources = 0 ;
    } else if (lParticleType == "Nuclei") {  //Changed by NN
      _fIsPhotonSources = 0 ;	
      int lMassNumber, lChargeNumber;
      double lAbundance;
      _fInitialChargeNumber.clear();
      _fInitialMassNumber.clear();
      _fAbundance.clear();
      TiXmlElement* lpNumberOfSpecies = lpXmlParticleType->FirstChildElement("Number_Of_Species") ;
      if(lpNumberOfSpecies) {
	lpNumberOfSpecies->Attribute("value", &_fNSpecies);
	TiXmlElement* lpSpecies = lpXmlParticleType->FirstChildElement("Species") ;
	for(int j=0; j< _fNSpecies; j++){
	  const char* lTst=lpSpecies->Attribute("MassNumber",&lMassNumber) ;
	  if(!lTst || lMassNumber <= 0 || lMassNumber > 56)throw TXmlErr("Error parsing MassNumber");
	  lTst=lpSpecies->Attribute("ChargeNumber", &lChargeNumber);
	  if(!lTst || lChargeNumber < 0 || lChargeNumber > 26)throw TXmlErr("Error parsing ChargeNumber");
	  lTst=lpSpecies->Attribute("Abundance", &lAbundance);
	  if(!lTst) lAbundance=1;
	  _fInitialMassNumber.push_back(lMassNumber);
	  _fInitialChargeNumber.push_back(lChargeNumber);
	  _fAbundance.push_back(lAbundance);
	  lpSpecies = lpSpecies->NextSiblingElement("Species") ;
	  if(!lpSpecies && j != _fNSpecies -1 ) throw TXmlErr("Error parsing nuclei species with the specified number in xml file" );

	}
      } else {
	_fNSpecies=0;
 	TiXmlElement* lpMassNumber = lpXmlParticleType->FirstChildElement("MassNumber") ;
 	if (!lpMassNumber) throw TXmlErr("No MassNumber for point source." );
 	lpMassNumber->Attribute("value",&lMassNumber) ;
 	_fInitialMassNumber.push_back(lMassNumber);
 	TiXmlElement* lpChargeNumber = lpXmlParticleType->FirstChildElement("ChargeNumber") ;
 	if (!lpChargeNumber) throw TXmlErr("No ChargeNumber for point source." );
 	lpChargeNumber->Attribute("value",&lChargeNumber) ;
 	_fInitialChargeNumber.push_back(lChargeNumber);
 	_fAbundance.push_back(1);
 	if(_fInitialChargeNumber.size()!=1 || _fInitialMassNumber.size() !=1 || _fAbundance.size() != 1 )
	TXmlErr("Error parsing single nuclei" );
      }
      double lAbundanceSum=0;
      for(int i=0; i < _fAbundance.size(); i++) lAbundanceSum+=_fAbundance[i];
      for(int i=0; i < _fAbundance.size(); i++) _fAbundance[i]/=lAbundanceSum;
      if(_fAbundance.size() > 1) for(int i=0; i < _fAbundance.size(); i++) _fAbundance[i]+=_fAbundance[i-1];

   } else throw TXmlErr("The injected particle type doesn't exist.") ;
  } else  throw TXmlErr("You haven't defined a particle type to be injected.");

  // 1 ) Number of sources
  TiXmlElement* lpXmlNumber = lpXmlSources->FirstChildElement("Number") ;
  if ( ! lpXmlNumber ) throw TXmlErr( "Error getting number of sources.");
  int lNb ;
  const char* lDummy = lpXmlNumber->Attribute("value",&lNb) ;
  if ( ! lDummy || lNb <= 0 )
    throw TXmlErr( "Error getting number of sources.");
  _fNumber = lNb ;
  
  // 2 ) Positions
  TiXmlElement *lpXmlCoordX, *lpXmlCoordY, *lpXmlCoordZ ;
  string lDum ;
  double lX,lY,lZ ;
  TVector3D lPos ;
  TiXmlElement* lpXmlDensity = lpXmlSources->FirstChildElement("Density") ;
  bool lIndividualEmax = 0 ; // set to 1 if Emax is set for individual sources
  if (lpXmlDensity) { // Density from which to get sources is specified

    TiXmlElement* lpPointSource = lpXmlSources->FirstChildElement("PointSource") ;
    if (lpPointSource) throw TXmlErr("Point sources are specified at the same time as a source density." );
    _fpSourceDensity = new TSourceDensity(aFileName) ;
    for (unsigned int j=0; j<_fNumber; j++)
      _fPositions.push_back(_fpSourceDensity->getSourcePosition()) ;
    
  } else { // point source list is specified
    
    _fpSourceDensity = NULL ;
    TiXmlElement* lpPointSource = lpXmlSources->FirstChildElement("PointSource") ;
    if (!lpPointSource) throw TXmlErr("No point sources, neither source density found." );
    TiXmlElement* lpXmlEmax = lpPointSource->FirstChildElement("Ecut_EeV") ;
    if(!lpXmlEmax){
      lpXmlEmax =lpPointSource->FirstChildElement("Rigidity_EeV") ;
      _fRigidityFlag=1;
    } else _fRigidityFlag=0;
    if (lpXmlEmax) lIndividualEmax = 1 ;
    for (unsigned int j=0; j<_fNumber; j++) {
      lpXmlCoordX = lpPointSource->FirstChildElement("CoordX_Mpc") ;
      if (!lpXmlCoordX) throw TXmlErr("No CoordX for point source.");
      lDum = lpXmlCoordX->Attribute("value",&lX) ;
      lPos.setX(lX*Mpc) ;	
      if (lUnivType == "One Dimension") {
	lPos.setY(0) ;
	lPos.setZ(0) ;
#ifdef DEBUG_OUTPUT 
	std::cout<< "Read discrete source (case 1d) # "<<j<<" at ("
		 <<lPos.x()<<","
		 <<lPos.y()<<","
		 <<lPos.z()<<")"
		 <<std::endl;
#endif
      } else {
	lpXmlCoordY = lpPointSource->FirstChildElement("CoordY_Mpc") ;
	if (!lpXmlCoordY) throw TXmlErr("No CoordY for point source." );
	lDum = lpXmlCoordY->Attribute("value",&lY) ;
	lpXmlCoordZ = lpPointSource->FirstChildElement("CoordZ_Mpc") ;
	if (!lpXmlCoordZ) throw TXmlErr("No CoordZ for point source." );
	lDum = lpXmlCoordZ->Attribute("value",&lZ) ;
	lPos.setY(lY*Mpc) ;
	lPos.setZ(lZ*Mpc) ;
#ifdef DEBUG_OUTPUT 
	std::cout<< "Read discrete source (case 3d) # "<<j<<" at ("
		 <<lPos.x()/Mpc<<","
		 <<lPos.y()/Mpc<<","
		 <<lPos.z()/Mpc<<")"
		 <<" / F.Y.I. unit conversion via Mpc="<<Mpc
		 <<std::endl;
#endif
      }
      // Get also the individual Emax if needed
      if (lIndividualEmax) {
	lpXmlEmax = lpPointSource->FirstChildElement("Ecut_EeV") ;
	if (!lpXmlEmax) {
	lpXmlEmax = lpPointSource->FirstChildElement("Rigidity_EeV") ;
	  _fRigidityFlag=1;
	} else _fRigidityFlag=0;
	if (!lpXmlEmax) throw TXmlErr("No individual Emax for all the sources") ;
	lDum = lpXmlEmax->Attribute("value",&_fEcut) ;
	_fEcut *= EeV ;
	_fEcutList.push_back(_fEcut); 
      }
      lpPointSource = lpPointSource->NextSiblingElement("PointSource") ;
      if ( (j!=_fNumber-1 && ! lpPointSource) || (j==_fNumber-1 && lpPointSource) ) 
	throw TXmlErr("Error parsing point sources with the specified number in xml file" );
      _fPositions.push_back(lPos) ;
    }

  }
  
  // 3 ) Spectrum
  TiXmlElement* lpSpectrum = lpXmlSources->FirstChildElement("Spectrum") ;
  if (!lpSpectrum) throw TXmlErr( "No source spectrum specified.");
  string lSpecType = lpSpectrum->Attribute("type") ;
  if (lSpecType == "Monochromatic") {
    _fSpectrumFlag = 1 ;
    if (!lIndividualEmax) {
      TiXmlElement* lpEnergy = lpSpectrum->FirstChildElement("Energy_EeV") ;
      if (!lpEnergy) {
	lpEnergy = lpSpectrum->FirstChildElement("Rigidity_EeV") ;
      if (!lpEnergy) throw TXmlErr("Neither Energy nor Rigidity in monochromatic spectrumspecified." );
      _fRigidityFlag=1;
      } else _fRigidityFlag=0;
      lDum = lpEnergy->Attribute("value",&_fEcut) ;
      _fEcut *= EeV ;
    }
    _fAlpha = 0 ;
    _fSigAlpha = 0 ;

  } else if (lSpecType == "Power Law") {

    _fSpectrumFlag = 2 ;
    TiXmlElement* lpAlpha = lpSpectrum->FirstChildElement("Alpha") ;
    if (!lpAlpha) throw TXmlErr("Spectral index unspecified." );
    lDum = lpAlpha->Attribute("value",&_fAlpha) ;
    TiXmlElement* lpSigAlpha = lpSpectrum->FirstChildElement("SigmaAlpha") ;
    if (!lpSigAlpha) _fSigAlpha = 0 ;
    else lDum = lpSigAlpha->Attribute("value",&_fSigAlpha) ;
    if (!lIndividualEmax) {
      TiXmlElement* lpEcut = lpSpectrum->FirstChildElement("Ecut_EeV") ;
       if (!lpEcut) {
      lpEcut = lpSpectrum->FirstChildElement("Rigidity_EeV") ;
      _fRigidityFlag=1;
      } else _fRigidityFlag=0;

      if (!lpEcut) throw TXmlErr("Spectral cut-off energy unspecified,");
      lDum = lpEcut->Attribute("value",&_fEcut) ;
      _fEcut *= EeV ;
      TiXmlElement* lpEmin = lpSpectrum->FirstChildElement("MinEnergy_EeV_Source") ;
//      if(!lpEmin){
//        lpEmin =  lpSpectrum->FirstChildElement("MinRigidity_EeV_Source") ;
//      }
      if(lpEmin){
      	lDum = lpEmin->Attribute("value",&_fEmin) ;
      	_fEmin*= EeV;
      }
    }

  } else throw TXmlErr("Error getting source spectrum type.");

  if (!lIndividualEmax) { 
    _fEcutList.resize(_fNumber);
    fill_n(_fEcutList.begin(), _fNumber, _fEcut);
    _fEminList.resize(_fNumber);
    fill_n(_fEminList.begin(), _fNumber, _fEmin);
  }
  double *lpFluct = new double[_fNumber];
  RandFlat::shootArray(_fNumber, lpFluct) ; // top-hat distribution
  for (unsigned int k=0; k<_fNumber; k++) 
    _fAlphaList.push_back(_fAlpha + _fSigAlpha*(2*lpFluct[k] - 1)) ;
  delete[] lpFluct ;

  // Consistency check
  if (!lIndividualEmax) {
    if (_fEcut <= _fEmin) throw TCrpErr("Minimum energy too large.") ;
  } else {
    for (unsigned int j=0; j<_fNumber; j++)
      if (_fEcutList.at(j) <= _fEmin) throw TCrpErr("Minimum energy too large.") ;
  }

}
