
/**
   @file    basicparam.cc
   @author  Tristan Beau, beau@in2p3.fr
   @brief   Implementation of the TBasicParam class. See the .h file
*/

#include "basicparam.h"
#include "TabulatedTALYSAveragedCrossSection.h"
#include "TabulatedTALYSMeanFreePath.h"

TBasicParam::TBasicParam(const char* aFileName) : TXmlParam(aFileName) {
  // Read number of particles to record
  _fN = XmlExtract().GetULong("TrajNumber") ;
  if (_fN <= 0 ) throw TXmlErr("Inappropriate Trajectory Number.") ;

  // Read minimal energy under which particles are abandonned
  _fEmin = XmlExtract().GetDouble("MinEnergy_EeV") ;
  _fEmin *= EeV ;

  // Read maximum time above which particles are abandonned
  _fTmax = XmlExtract().GetDouble("MaxTime_Mpc") ;
  _fTmax *= Mpc * inv_c_light ;


  // Possibility to set the seed for random numbers
  bool lCheck = XmlExtract().CheckElement("RandomSeed") ;
  if (lCheck) {
    _fHepRdmSeed = (long)(XmlExtract().GetInt("RandomSeed")) ;
    RandFlat::setTheSeed(_fHepRdmSeed) ;
  } else _fHepRdmSeed = RandFlat::getTheSeed() ;

  // Output type
  TiXmlElement* lpXmlOutput = XmlExtract().GetElement("Output") ;
  string lRecordType ;
  try {
    lRecordType = GetAttribute(lpXmlOutput,"type") ;
  }
  catch ( TXmlErr ) {
    throw TXmlErr ("Output type undefined.");
  }
  if (lRecordType == "None" ||
      lRecordType == "Full Trajectories" ||
      lRecordType == "Events") _fRecordMode = lRecordType ;
  else throw TXmlErr("Bad determination of Recording mode" );

  // Output file name and format
  if (_fRecordMode != "None") {
    TiXmlElement* lpXmlFile = lpXmlOutput->FirstChildElement("File") ;
    if ( ! lpXmlFile ) throw TXmlErr("Output file not defined");
    TiXmlNode* lpXmlFileName = lpXmlFile->FirstChild() ;
    if ( (! lpXmlFileName) || (! lpXmlFileName->ToText()) ) 
      throw TXmlErr("Bad output file name") ;
    _fOutputFileName = lpXmlFileName->Value() ;


    // Overwrite mode or not
    bool lForce = false;  
    string lTemp ;
    try {
      lTemp = GetAttribute(lpXmlFile,"option") ;
      if ( lTemp == "force" )
	lForce = true;
      else
	throw TCrpErr ( string("Unknown option used for output file : ") + lTemp );
    }
    catch ( TXmlErr ) {
      ; // option undefined. Use default, i.e. lForce = false
    }
    
    
    string lOutputFormat = lpXmlFile->Attribute("type") ;
    if (lOutputFormat == "ASCII") {
      _fpOutputData = new TAsciiOutput(aFileName, _fOutputFileName, lForce) ;
    } else if (lOutputFormat == "FITS") {
      _fpOutputData = new TFitsOutput(aFileName, _fOutputFileName, lForce) ;
    } else if (lOutputFormat == "ROOT") {
      _fpOutputData = new TRootOutput(aFileName, _fOutputFileName, lForce) ;
    } else if (lOutputFormat == "HDF") {
      throw TCrpErr("HDF output format not yet implemented");
    } else {
      throw TXmlErr("Bad output file format");
    }

  } else {
    _fpOutputData = new TNoOutput ;
  }

} 

TBasicParam::~TBasicParam() {
  delete _fpOutputData ;
}
