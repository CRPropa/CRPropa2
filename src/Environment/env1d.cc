
/** (beta)
    @file    env1d.cc
    @author  Eric Armengaud, armengau@in2p3.fr
    @brief   Implementation of the TEnv1D class. See the .h file
*/

#include "env1d.h"
#include "CMBLikeScaling.h"
#include "NoIRBScaling.h"
#include "UserDefinedScaling.h"
#include "Kneiske2004_BestFit.h"

TEnv1D::TEnv1D(const char* aFileName) : TXmlParam(aFileName) {
  SetType(UNIVERSE_ENV1D);


  //Read settings from xml configuration
  //TXmlParam Test(XMLFileName.c_str());
  _fModel = 0;
  bool lCheckIRBzModel = XmlExtract().CheckElement("IRB_MFPScalingModel") ;
  if ( ! lCheckIRBzModel){
    std::cout<<"No IRB evolution scaling model chosen via xml settings. Default \"cmblike\"-scaling will be used."<<std::endl;
    _fModel = new CMBLikeScaling(aFileName);
  }
  else {
    TiXmlElement* lpXmlIRBzModel = XmlExtract().GetElement("IRB_MFPScalingModel") ;
    std::string lIRBzModelType = lpXmlIRBzModel->Attribute("type") ;
    std::cout<<"IRB_MFPScalingModel: "<<lIRBzModelType<<"is selected."<<std::endl;
    if (lIRBzModelType == "CMBLike") {
      std::cout<<"IRB_MFPScalingModel: "<<lIRBzModelType<<"is selected."<<std::endl;
      _fModel = new CMBLikeScaling(aFileName);
    }
    else if (lIRBzModelType == "Kneiske2004_BestFit") {
      std::cout<<"IRB_MFPScalingModel: "<<lIRBzModelType<<" is selected."<<std::endl;
      _fModel = new Kneiske2004_BestFit(aFileName);
    }
    else if (lIRBzModelType == "UserDefinedScaling") {
          std::cout<<"IRB_MFPScalingModel: "<<lIRBzModelType<<" is selected."<<std::endl;
          _fModel = new UserDefinedScaling(aFileName);
        }
    else if (lIRBzModelType == "None") {
      std::cout<<"IRB_MFPScalingModel: "<<lIRBzModelType<<" is selected."<<std::endl;
      _fModel = new NoIRBScaling(aFileName);
    }
    else throw TCrpErr("<IRB_MFPScalingModel type=XYZ> selected but type not specified.") ;
  }

  // Magnetic field
  _fBflag = 0; // even if there is a magnetic field for photons,
  // we don't want any deflection!
  bool lCheckMag = XmlExtract().CheckElement("MagneticField") ;
  if ( !lCheckMag ) _fpMagField = new TNullMagField() ;
  else {
    TiXmlElement* lpXmlField = XmlExtract().GetElement("MagneticField") ;
    string lFieldType = lpXmlField->Attribute("type") ;
    if (lFieldType == "1D") _fpMagField = new TField1D(aFileName) ;
    else if(lFieldType == "Null") _fpMagField = new TNullMagField();
    else throw TCrpErr("Incorrect field type for 1D Environment") ;
  }

  // Interactions
  TiXmlElement* lpXmlInteract = XmlExtract().GetElement("Interactions") ;
  string lInteractionType = lpXmlInteract->Attribute("type") ;
  if (lInteractionType == "F77-proton") {
    _fpInteractionData = new TBasicPInteractions(aFileName) ;
  } else if (lInteractionType == "None") {
    _fpInteractionData = new TNullInteractions() ;
  } else if (lInteractionType == "Sophia") {
    _fpInteractionData = new TSophiaInteractions(aFileName,_fModel) ;
  } else if (lInteractionType == "Photon") {
    _fpInteractionData = new TPhotonInteractions(aFileName) ;
  } else throw TXmlErr("Error getting the interaction type.");

  // Observer
  bool lCheckObs = XmlExtract().CheckElement("Observers") ;
  if (lCheckObs == 1) 
    throw TXmlErr("No observer should be present in 1D");
  TiXmlElement* lpXmlOutput = XmlExtract().GetElement("Output") ;
  string lRecordType = lpXmlOutput->Attribute("type") ;
  if (lRecordType == "None" || lRecordType == "Full Trajectories") {
    _fpObservers = new TNoObserver() ;
  } else if (lRecordType == "Events") {
    _fpObservers = new TPointObserver(aFileName) ;
  } else throw TXmlErr("Bad determination of output type" );  

  // Sources
  TiXmlElement* lpXmlSources = XmlExtract().GetElement("Sources") ;
  string lSourceType = lpXmlSources->Attribute("type") ;
  if (lSourceType == "Discrete") {
    _fpSources = new TDiscreteSources(aFileName) ; 
  } else if (lSourceType == "Continuous") {
    _fpSources = new TContinuousSources(aFileName) ;
  } else throw TXmlErr("Unknown source type");

  // Backgrounds
  bool lCheckIR = XmlExtract().CheckElement("InfraredBackground") ;
  string lIRType ;
  if (lCheckIR) {
    TiXmlElement* lpXmlIR = XmlExtract().GetElement("InfraredBackground") ;
    lIRType = lpXmlIR->Attribute("type") ;
  } else lIRType = "Uniform" ;
  if (lIRType == "Uniform") {
    _fpInfrared = new TUniformInfrared(aFileName) ;
  } else if (lIRType == "Variable") { // LM
    _fpInfrared = new TVariableInfrared(aFileName, true);
  } else throw TXmlErr("Unknown IR background type") ;
  // Not yet implemented grid IR :-)

  bool lCheckRadio = XmlExtract().CheckElement("RadioBackground") ;
  // Accepted values are "Obs", "Med", "High", "Null"
  if (lCheckRadio) {
    TiXmlElement* lpXmlRadio = XmlExtract().GetElement("RadioBackground") ;
    _fRadioBackground = lpXmlRadio->Attribute("type") ;
  } else _fRadioBackground = "Obs" ;
  if (_fRadioBackground != "High" && _fRadioBackground != "Med" && _fRadioBackground != "Obs" && _fRadioBackground != "Null")
    throw TXmlErr("Unknown Radio background type") ;

  // Gas
  bool lCheckGas = XmlExtract().CheckElement("Gas") ;
  string lGasType = "" ;
  if (lCheckGas) {
    TiXmlElement* lpXmlGas = XmlExtract().GetElement("Gas") ;
    lGasType = lpXmlGas->Attribute("type") ;
  } 
  if (lGasType == "SpherSymm") _fpGas = new TClusterGas(aFileName) ;
  else {
    _fpGas = NULL;
    //cerr << "Requested 1D Gas type not yet implemented" << endl;
  }

  // Distance-Redshift array
  bool lCheckOmegaM = XmlExtract().CheckElement("OmegaM") ;
  if (lCheckOmegaM) {
    _fOmegaM = XmlExtract().GetDouble("OmegaM") ;
  } else _fOmegaM = DEFAULT_OMEGA_M ;
  bool lCheckOmegaLambda = XmlExtract().CheckElement("OmegaLambda") ;
  if (lCheckOmegaLambda) {
    _fOmegaLambda = XmlExtract().GetDouble("OmegaLambda") ;
  } else _fOmegaLambda = DEFAULT_OMEGA_LAMBDA ;
  bool lCheckH0 = XmlExtract().CheckElement("H0_km_s_Mpc") ;
  if (lCheckH0) {
    _fH0 = XmlExtract().GetDouble("H0_km_s_Mpc") ;
  } else _fH0 = DEFAULT_H_0_KM_S_MPC ;
  _fH0 *= kilometer/(second*Mpc) ;

  long lBinsArrays = 1000 ; // tout ca en dur pour l'instant...
  double lZmin=0.0001 ;
  double lZmax = 100 ;
  _fRedshiftArray.push_back(0) ;
  for (long i=0; i < lBinsArrays-1; i++) {
    double kk=i/(lBinsArrays-2.) ; 
    _fRedshiftArray.push_back(lZmin*pow(lZmax/lZmin,kk)) ;
  }

  vector<double> lEvolutionFactor ;
  lEvolutionFactor.push_back(1) ;
  _fDistanceArray.push_back(0) ;
  for (unsigned long i=1; i < _fRedshiftArray.size(); i++ ) {
    double lZ = _fRedshiftArray[i] ;
    lEvolutionFactor.push_back(1./((1.+lZ)*sqrt(_fOmegaLambda+_fOmegaM*pow(1.+lZ,3)))) ;
    _fDistanceArray.push_back(_fDistanceArray[i-1]+(lZ-_fRedshiftArray[i-1])*0.5*(lEvolutionFactor[i-1]+lEvolutionFactor[i])*c_light/_fH0) ;
  }

  // Xmax is determined from the sources
  if (lSourceType == "Discrete") {
    _fXmax = 0 ;
    for (unsigned int i=0; i< _fpSources->Nb(); i++) {
      if ( _fXmax < (_fpSources->Positions(i)).x() ) _fXmax = (_fpSources->Positions(i)).x() ;
    }
  } else _fXmax = _fpSources->Density()->Xmax() ;
  if ( _fXmax  > _fDistanceArray.back() ) throw TCrpErr("Sources more distant than maximum redshift.\nMaximum distance is 4 Gpc.") ;
  _fXmax *= 1.01 ;

  // Tests
  if ( lSourceType == "Discrete" && lCheckMag ) {
    for (unsigned int i=0; i< _fpSources->Nb(); i++) {
      if ( _fpMagField->Xmax() < (_fpSources->Positions(i)).x() )
	throw TCrpErr("B field Xmax not large enough for the sources.") ;
    }
  }

  // more .. (todo...) check same for continuous sources...

}

TEnv1D::~TEnv1D() {

  delete _fpMagField ;
  delete _fpInteractionData ;
  delete _fpObservers ;
  delete _fpSources ;
  delete _fpInfrared ;
  delete _fpGas;
  delete _fModel;
}
