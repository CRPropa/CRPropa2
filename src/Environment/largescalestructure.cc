
/**(beta)
   @file    largescalestructure.cc
   @author  Eric Armengaud, armengau@in2p3.fr
   @brief   Implementation of the TLargeScaleStructure class. See the .h file
*/

#include "largescalestructure.h"
#include "CMBLikeScaling.h"
#include "NoIRBScaling.h"
#include "Kneiske2004_BestFit.h"
#include "UserDefinedScaling.h"

TLargeScaleStructure::TLargeScaleStructure(const char* aFileName) : TXmlParam(aFileName){
	SetType(UNIVERSE_LARGESCALE) ;

  _fModel = 0;
  bool lCheckIRBzModel = XmlExtract().CheckElement("IRB_MFPScalingModel") ;
  if ( ! lCheckIRBzModel){
    std::cout<<"No IRB evolution scaling model chosen via xml settings. Default \"cmblike\"-scaling will be used."<<std::endl;
    _fModel = new CMBLikeScaling(aFileName);
  }
  else {
    TiXmlElement* lpXmlIRBzModel = XmlExtract().GetElement("IRB_MFPScalingModel") ;
    std::string lIRBzModelType = lpXmlIRBzModel->Attribute("type") ;
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
    
	string lDum ;
  double lBx, lBy, lBz ;

  // Bounds for the simulation box
  TiXmlElement* lpEnv = XmlExtract().GetElement("Environment") ;
  TiXmlElement* lpXmin = lpEnv->FirstChildElement("Xmin_Mpc") ;
  TiXmlElement* lpXmax = lpEnv->FirstChildElement("Xmax_Mpc") ;
  TiXmlElement* lpYmin = lpEnv->FirstChildElement("Ymin_Mpc") ;
  TiXmlElement* lpYmax = lpEnv->FirstChildElement("Ymax_Mpc") ;
  TiXmlElement* lpZmin = lpEnv->FirstChildElement("Zmin_Mpc") ;
  TiXmlElement* lpZmax = lpEnv->FirstChildElement("Zmax_Mpc") ;
  if (lpXmin) {
    if (!lpXmax || !lpYmin || !lpYmax || !lpZmin || !lpZmax)
      throw TXmlErr("Environment boundary missing.") ;
    lDum = lpXmin->Attribute("value",&_fXmin) ;
    _fXmin *= Mpc ;
    lDum = lpXmax->Attribute("value",&_fXmax) ;
    _fXmax *= Mpc ;
    lDum = lpYmin->Attribute("value",&_fYmin) ;
    _fYmin *= Mpc ;
    lDum = lpYmax->Attribute("value",&_fYmax) ;
    _fYmax *= Mpc ;
    lDum = lpZmin->Attribute("value",&_fZmin) ;
    _fZmin *= Mpc ;
    lDum = lpZmax->Attribute("value",&_fZmax) ;
    _fZmax *= Mpc ;
  }

  // Magnetic field
  TiXmlElement* lpXmlField = XmlExtract().GetElement("MagneticField") ;
  string lFieldType = lpXmlField->Attribute("type") ;
  if (lFieldType == "Null") {

    if (!lpXmin) throw TXmlErr("Environment boundaries not specified.") ;
    _fBflag = 0 ;
    _fpMagField = new TNullMagField() ;

  } else if (lFieldType == "Uniform") {

    if (!lpXmin) throw TCrpErr("Environment boundaries not specified.") ;
    _fBflag = 1 ;
    TiXmlElement* lpXmlBx = lpXmlField->FirstChildElement("Bx_nG") ;
    TiXmlElement* lpXmlBy = lpXmlField->FirstChildElement("By_nG") ;
    TiXmlElement* lpXmlBz = lpXmlField->FirstChildElement("Bz_nG") ;
    if (!lpXmlBx || !lpXmlBy || !lpXmlBz) {
      throw TCrpErr("Magnetic field value not given in xml file"); 
    }
    lDum = lpXmlBx->Attribute("value",&lBx) ;
    lDum = lpXmlBy->Attribute("value",&lBy) ;
    lDum = lpXmlBz->Attribute("value",&lBz) ;
    TVector3D lField(lBx*nG,lBy*nG,lBz*nG) ;
    _fpMagField = new TUniformMagField(lField) ;

  } else if (lFieldType == "LSS-Grid") {
    _fBflag = 1 ;
    _fpMagField = new TLSSField(aFileName) ;
    if (lpXmin || lpXmax || lpYmin || lpYmax || lpZmin || lpZmax)
      throw TXmlErr("X/Y/Zmin/max should not be defined with magnetic grid");
    _fXmin = _fpMagField->Origin().x() ;
    _fXmax = _fXmin + _fpMagField->Stepsize()*_fpMagField->Nx() ;
    _fYmin = _fpMagField->Origin().y() ;
    _fYmax = _fYmin + _fpMagField->Stepsize()*_fpMagField->Ny() ;
    _fZmin = _fpMagField->Origin().z() ;
    _fZmax = _fZmin + _fpMagField->Stepsize()*_fpMagField->Nz() ;

  } else if (lFieldType == "Kolmogoroff") {
    // Kolmogoroff field
    _fBflag = 1 ;
    _fpMagField = new TKolmogoroffMagField(aFileName) ;
    if (lpXmin || lpXmax || lpYmin || lpYmax || lpZmin || lpZmax) 
      throw TXmlErr("X/Y/Zmin/max should not be defined with magnetic grid") ;
    _fXmin = _fpMagField->Origin().x() ;
    _fXmax = _fXmin + _fpMagField->Stepsize()*_fpMagField->Nx() ;
    _fYmin = _fpMagField->Origin().y() ;
    _fYmax = _fYmin + _fpMagField->Stepsize()*_fpMagField->Ny() ;
    _fZmin = _fpMagField->Origin().z() ;
    _fZmax = _fZmin + _fpMagField->Stepsize()*_fpMagField->Nz() ;

    //     Output of the magnetic field if requested
    // This should be moved to gridmagnetic field
    // Not yet done, because the only function called is the constructor
    // And then the field is not created. 

    TiXmlElement* outFile = lpXmlField->FirstChildElement("Outputfile");

    if(outFile){//Output of the magnetic field. 
      string outFOption;
      string outFName(outFile->Attribute("name"));
      try{//GetAttribute raises an exception, if "option" is not present.
	outFOption=GetAttribute(outFile, "option");
	if(outFOption=="force"){
	  outFName="!"+outFName; // Inherited from fits_create_file(), which assumes that a filename starting with '!'
	  //means overwriting of an existing file, cropping the '!'. see cfitsio documentation. 
	} else throw TCrpErr( string("Unknown option used for output file : ") + outFOption); 
      } catch (TXmlErr){
	outFOption="none";
	// outFOption not set. 
      };

      if(_fpMagField->writeMagField(outFName)){//
	std::cout<< "Output magnetic field " << outFName << "written." << std::endl;
      } else {
	std::cout<< "No magnetic field output written." << std::endl;
	// This should never be triggered because writeMagField should throw an error. 
	//( It would be triggered if the virtual function TMagField::writeMagField is called instead of TGridField::writeMagField.)
      }
    }

  } else throw TXmlErr("No correct magnetic field specified");

  // Gas
  bool lCheckGas = XmlExtract().CheckElement("Gas") ;
  string lGasType ;
  if (lCheckGas) {
    TiXmlElement* lpXmlGas = XmlExtract().GetElement("Gas") ;
    lGasType = lpXmlGas->Attribute("type") ;
    if (lGasType == "SpherSymm") _fpGas = new TClusterGas(aFileName) ;
    else if (lGasType == "LSS-Grid") {
      _fpGas = new TLSSGas(aFileName) ;
      
      if (lpXmin || lpXmax || lpYmin || lpYmax || lpZmin || lpZmax)
	throw TXmlErr("X/Y/Zmin/max should not be defined with gas grid");
      _fXmin = _fpGas->Origin().x() ;
      _fXmax = _fXmin + _fpGas->Stepsize()*_fpGas->Nx() ;
      _fYmin = _fpGas->Origin().y() ;
      _fYmax = _fYmin + _fpGas->Stepsize()*_fpGas->Ny() ;
      _fZmin = _fpGas->Origin().z() ;
      _fZmax = _fZmin + _fpGas->Stepsize()*_fpGas->Nz() ;
      
    } 
  } else {
    _fpGas = NULL; 
    //cerr << "Requested 3D Gas type not yet implemented" << endl;
  }
  // Observers
  TiXmlElement* lpXmlOutput = XmlExtract().GetElement("Output") ;
  string lRecordType = lpXmlOutput->Attribute("type") ;
  if (lRecordType == "None" || lRecordType == "Full Trajectories") {
    _fpObservers = new TNoObserver() ;
    //    if (lpXmlObs) throw TXmlErr("There should not be any Observer.") ;
  } else if (lRecordType == "Events") {
    TiXmlElement* lpXmlObs = XmlExtract().GetElement("Observers") ;
    if (!lpXmlObs) throw TXmlErr("There should be an Observer.") ;
    lDum = lpXmlObs->Attribute("type") ;
    if (lDum == "Spheres around Observers")
      _fpObservers = new TSmallSphereObservers(aFileName) ;
    else if (lDum == "Spheres around Source")
      _fpObservers = new TLargeSphereObservers(aFileName) ;
    else throw TXmlErr("Observer type unspecified.") ;
    for (int i=0; i<_fpObservers->Nb(); i++) {
      TVector3D lPos = _fpObservers->Positions(i) ;
      if (lPos.x() > _fXmax || lPos.x() < _fXmin ||
	  lPos.y() > _fYmax || lPos.y() < _fYmin ||
	  lPos.z() > _fZmax || lPos.z() < _fZmin)
	throw TCrpErr("Observer positions not inside the simulation box.") ;
    } 
    
  } else throw TXmlErr("Bad determination of output type" ); 

  // Sources
  TiXmlElement* lpXmlSources = XmlExtract().GetElement("Sources") ;
  string lSourceType = lpXmlSources->Attribute("type") ;
  if (lSourceType == "Discrete") {
    _fpSources = new TDiscreteSources(aFileName) ;
    for (unsigned int i=0; i<_fpSources->Nb(); i++) {
      TVector3D lPos = _fpSources->Positions(i) ;
      if (lPos.x() > _fXmax || lPos.x() < _fXmin ||
	  lPos.y() > _fYmax || lPos.y() < _fYmin ||
	  lPos.z() > _fZmax || lPos.z() < _fZmin)
	throw TCrpErr("Source positions not inside the simulation box.") ;
    }
  } else if (lSourceType == "Continuous") {
    _fpSources = new TContinuousSources(aFileName) ;
  } else throw TCrpErr("Error getting the source type in xml file");

  // Interactions
  TiXmlElement* lpXmlInteract = XmlExtract().GetElement("Interactions") ;
  string lInteractionType = lpXmlInteract->Attribute("type") ;
  if (lInteractionType == "F77-proton") {
    _fpInteractionData = new TBasicPInteractions(aFileName) ;
  } else if (lInteractionType == "None") {
    _fpInteractionData = new TNullInteractions() ;
  } else if (lInteractionType == "Sophia") {
    _fpInteractionData = new TSophiaInteractions(aFileName, _fModel) ;
  } else if (lInteractionType == "Photon") {
    _fpInteractionData = new TPhotonInteractions(aFileName) ;
  } else throw TCrpErr("Error getting the interaction type in xml file"); 

  // Integrator for deflection equations
  TiXmlElement* lpXmlIntegrator = XmlExtract().GetElement("Integrator") ;
  _fIntegrator = lpXmlIntegrator->Attribute("type") ;
  if (_fIntegrator == "Cash-Karp RK") {
    TiXmlElement* lpXmlEps = lpXmlIntegrator->FirstChildElement("Epsilon") ;
    if (!lpXmlEps) {
      throw TCrpErr("Epsilon unspecified for RKCK integrator in xml file");
    }
    lDum = lpXmlEps->Attribute("value",&_fIntegratorEps) ;
  } else if (_fIntegrator == "Bulirsch-Stoer") {
    throw TCrpErr("Bulirsch-Stoer integrator not fully implemented.") ;
  } else throw TXmlErr("Incorrect integrator type.");

  _fIntegratorMinTimeStep = min(_fpInteractionData->InteractionTimeStep(),
				MINSTEP_FACTOR*_fpMagField->Stepsize()/c_light) ;
  // Possibility to fix by hand the minimum step for integrator
  // By default it will be defined by MINSTEP_FACTOR
  TiXmlElement* lpMinstep = lpXmlIntegrator->FirstChildElement("MinStep_Mpc") ;
  if (lpMinstep) {
    double lToto ;
    lDum = lpMinstep->Attribute("value",&lToto) ;
    _fIntegratorMinTimeStep = min(_fpInteractionData->InteractionTimeStep(),
				  lToto*Mpc/c_light) ;
  }
  // will not be 0 in the case of a grid field...

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
    _fpInfrared = new TVariableInfrared(aFileName,false);
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

  // Tests

  if (_fpObservers->Type() == OBSERVER_LARGESPHERE ) {
    if (_fpSources->Type() != SOURCE_DISCRETE || _fpSources->Nb() != 1) 
      throw TXmlErr("Only one discrete source with large spherical observers.") ;
    for (int j=0; j<_fpObservers->Nb(); j++) {
      if ( (_fpObservers->Positions)(j) != (_fpSources->Positions)(0) )
	throw TXmlErr("Spheres around source must have same position as source.") ; 
    }
    if ( _fpObservers->Nb() > 1 && (_fpInteractionData->SecPairProdPhotonFlag() ||
				     _fpInteractionData->SecPhotonFlag() || 
				     _fpInteractionData->SecNuFlag()))
      throw TXmlErr("Only one spherical observer is allowed for detection of secondaries.") ;
  }
  
  if (_fpSources->Density() != NULL) {
    TSourceDensity *lDens = _fpSources->Density() ;
    if (abs(lDens->Xmin()-_fXmin) > GRID_MATCH_MPC * Mpc || 
	abs(lDens->Xmax()-_fXmax) > GRID_MATCH_MPC * Mpc || 
    	abs(lDens->Ymin()-_fYmin) > GRID_MATCH_MPC * Mpc || 
	abs(lDens->Ymax()-_fYmax) > GRID_MATCH_MPC * Mpc || 
    	abs(lDens->Zmin()-_fZmin) > GRID_MATCH_MPC * Mpc || 
	abs(lDens->Zmax()-_fZmax) > GRID_MATCH_MPC * Mpc)
      throw TXmlErr("Source density grid do not match the environment boundaries") ;
  }

}

TLargeScaleStructure::~TLargeScaleStructure() {

  delete _fpMagField ;
  delete _fpObservers ;
  delete _fpSources ;
  delete _fpInteractionData ;
  delete _fpInfrared ;
  delete _fModel;

}
