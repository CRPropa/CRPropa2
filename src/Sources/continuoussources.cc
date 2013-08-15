/**
   @file   continuoussources.cc
   @author Eric Armengaud, armengau@in2p3.fr
   @brief  Implementation of the TContinuousSources class. See the .h file
*/

#include "continuoussources.h"

TContinuousSources::TContinuousSources(const char* aFileName) : TXmlParam(aFileName) {
  SetType(SOURCE_CONTINUOUS);

  //TODO: Parts of this XML parsing is done in TDiscreteSources as well?!? Here one can save some lines of code(!) : Create a class e.g.: TContinuousAndDiscreteSources with a suitable routine an inherit from there! 
  string ldum;
  TiXmlElement* lpXmlSources = XmlExtract().GetElement("Sources") ;
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

  // 1 ) Density
  string lDum ;
  TiXmlElement* lpXmlDensity = lpXmlSources->FirstChildElement("Density") ;
  if (!lpXmlDensity) throw TXmlErr("No source density specified for continuous sources.") ;
  _fpSourceDensity = new TSourceDensity(aFileName) ;

  // 2 ) Spectrum
  TiXmlElement* lpSpectrum = lpXmlSources->FirstChildElement("Spectrum") ;
  if (!lpSpectrum) throw TXmlErr( "Error : no spectrum specified in xml file");
  string lSpecType = lpSpectrum->Attribute("type") ;
  if (lSpecType == "Monochromatic") {

    _fSpectrumFlag = 1 ;
    TiXmlElement* lpEnergy = lpSpectrum->FirstChildElement("Energy_EeV") ;
      _fRigidityFlag=0;
    if(!lpEnergy){ lpEnergy = lpSpectrum->FirstChildElement("Rigidity_EeV") ;
      _fRigidityFlag=1;
    }

    if (!lpEnergy)
      throw TXmlErr("Neither Energy nor Rigidity in monochromatic spectrum specified." );
    lDum = lpEnergy->Attribute("value",&_fEcut) ;
#ifdef DEBUG_OUTPUT
    std::cout<<"TEST Energy from XML : Energy_EeV: "<<_fEcut
	     <<"\n\t\tGeV: "<<GeV
	     <<"\n\t\tEeV: "<<EeV
	     <<"\n\t\texaelectronvolt: "<<exaelectronvolt
	     <<"\n\t\tEnergy: "<<_fEcut*EeV
	     <<std::endl;
#endif
    _fEcut *= EeV ;

  } else if (lSpecType == "Power Law") {

    _fSpectrumFlag = 2 ;
    TiXmlElement* lpAlpha = lpSpectrum->FirstChildElement("Alpha") ;
    if (!lpAlpha) throw TXmlErr("Spectral index unspecified." );
    lDum = lpAlpha->Attribute("value",&_fAlpha) ;
    TiXmlElement* lpEcut = lpSpectrum->FirstChildElement("Ecut_EeV") ;
    _fRigidityFlag=0;
    if(!lpEcut){
      lpEcut =  lpSpectrum->FirstChildElement("Rigidity_EeV") ;
      _fRigidityFlag=1;
    }
    if (!lpEcut) throw TXmlErr("Spectral cut-off energy unspecified,");
    lDum = lpEcut->Attribute("value",&_fEcut) ;
    _fEcut *= EeV ;

    TiXmlElement* lpEmin = lpSpectrum->FirstChildElement("MinEnergy_EeV_Source") ;
//	if(!lpEmin){
//	  lpEmin =  lpSpectrum->FirstChildElement("MinRigidity_EeV_Source") ;
//	}
	if(lpEmin){
		lDum = lpEmin->Attribute("value",&_fEmin) ;
		_fEmin*= EeV;
	}
  } else throw TXmlErr("Error getting spectrum type.");

  // Consistency check
  if (_fEcut <= _fEmin) throw TCrpErr("Minimum energy too large.") ;

}

TContinuousSources::~TContinuousSources() {
  delete _fpSourceDensity ;
}
