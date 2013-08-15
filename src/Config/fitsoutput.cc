

/**
   @file   fitsoutput.cc
   @author Eric Armengaud, armengau@in2p3.fr
   @brief  Implementation of the TFitsOutput class. See the .h file
*/

#include "fitsoutput.h"

TFitsOutput::TFitsOutput(const char* aFileName, string aOutputFile, bool aForceOverwrite) : 
TXmlParam(aFileName) {
	SetType(OUTPUT_FITS);
  string lUnivType = (XmlExtract().GetElement("Environment"))->Attribute("type");
  string lRecordType = (XmlExtract().GetElement("Output"))->Attribute("type") ;
  TiXmlElement* lpXmlInter = XmlExtract().GetElement("Interactions") ;
  TiXmlElement* lpXmlCheckPhoton = lpXmlInter->FirstChildElement("SecondaryPhotons") ;
  TiXmlElement* lpXmlCheckPPPhoton = lpXmlInter->FirstChildElement("SecondaryPairProdPhotons") ;
  TiXmlElement* lpXmlCheckNu = lpXmlInter->FirstChildElement("SecondaryNeutrinos") ;
  string lInteractionType = lpXmlInter->Attribute("type") ;

  _fStatus = 0 ;
  if ( aForceOverwrite ) {
    string lRmCommand = "\\rm -f " + aOutputFile ;
    system ( lRmCommand.c_str() ) ;
  } 
  else {
    // test if file exists. If it does, throw exception because cannot recreate a fits file
    ifstream lFs(aOutputFile.c_str());
    if ( lFs ) throw TCrpErr("Fits output file exists. Cannot overwrite by default. Use option=\"force\" to overwrite it." );
    lFs.close();
  }
  
  fits_create_file(&_fpFitsStream, aOutputFile.c_str(), &_fStatus) ;
  if (_fStatus) throw TFitsErr(_fStatus) ;
  
  if (lUnivType == "One Dimension" && lRecordType == "Events") {
    _fNbFields = 6 ;
    char extname[] = "EVENT_TABLE" ;
    char *tableType[] = {(char*)"Particle_Type",
			 (char*)"Initial_Type", 
			 (char*)"Initial_Position",
			 (char*)"Initial_Energy",
			 (char*)"Time",
			 (char*)"Energy"} ;
    char *tableForm[] = {(char*)"1I",(char*)"1I",(char*)"1D",(char*)"1D",(char*)"1D",(char*)"1D"} ;
    char *tableUnit[] = {(char*)"None",(char*)"None",(char*)"Mpc",(char*)"EeV",(char*)"Mpc",(char*)"EeV"} ;
    fits_create_tbl(_fpFitsStream, BINARY_TBL, 0, _fNbFields, tableType, tableForm, tableUnit, extname, &_fStatus) ;
  } else if (lUnivType == "One Dimension" && lRecordType == "Full Trajectories") {
    _fNbFields = 5 ;
    char extname[] = "TRAJECTORY_TABLE" ;
    char *tableType[] = {(char*)"Particle_Type",
			 (char*)"Initial_Type",
			 (char*)"Time",
			 (char*)"Position",
			 (char*)"Energy"} ;
    char *tableForm[] = {(char*)"1I",(char*)"1I",(char*)"1D",(char*)"1D",(char*)"1D"} ;
    char *tableUnit[] = {(char*)"None",(char*)"None",(char*)"Mpc",(char*)"Mpc",(char*)"EeV"} ;
    fits_create_tbl(_fpFitsStream, BINARY_TBL, 0, _fNbFields, tableType, tableForm, tableUnit, extname, &_fStatus) ;
  } else if (lUnivType != "One Dimension" && lRecordType == "Events") {
    _fNbFields = 15 ;
    char extname[] = "EVENT_TABLE" ;
    char *tableType[] = {(char*)"Particle_Type",
			 (char*)"Initial_Type",
			 (char*)"Initial_Position_X",
			 (char*)"Initial_Position_Y",
			 (char*)"Initial_Position_Z",
			 (char*)"Initial_Energy",
			 (char*)"Initial_Momentum_Theta",
			 (char*)"Initial_Momentum_Phi",
			 (char*)"Time",
			 (char*)"Position_X",
			 (char*)"Position_Y",
			 (char*)"Position_Z",
			 (char*)"Energy",
			 (char*)"Momentum_Theta",
			 (char*)"Momentum_Phi"} ;
    char *tableForm[] = {(char*)"1I",(char*)"1I",(char*)"1D",(char*)"1D",(char*)"1D",(char*)"1D",(char*)"1D",(char*)"1D",(char*)"1D",(char*)"1D",(char*)"1D",(char*)"1D",(char*)"1D",(char*)"1D",(char*)"1D"} ;
    char *tableUnit[] = {(char*)"None",(char*)"None",(char*)"Mpc",(char*)"Mpc",(char*)"Mpc",(char*)"EeV",(char*)"Rad",(char*)"Rad",(char*)"Myr",(char*)"Mpc",(char*)"Mpc",(char*)"Mpc",(char*)"EeV",(char*)"Rad",(char*)"Rad"} ;
    fits_create_tbl(_fpFitsStream, BINARY_TBL, 0, _fNbFields, tableType, tableForm, tableUnit, extname, &_fStatus) ;
  } else if (lUnivType != "One Dimension" && lRecordType == "Full Trajectories") {
    _fNbFields = 10 ;
    char extname[] = "TRAJECTORY_TABLE" ;
    char *tableType[] = {(char*)"Particle_Type",
			 (char*)"Initial_Type",
			 (char*)"Time",
			 (char*)"Position_X",
			 (char*)"Position_Y",
			 (char*)"Position_Z",
			 (char*)"Momentum_X",
			 (char*)"Momentum_Y",
			 (char*)"Momentum_Z",
			 (char*)"Energy"} ;
    char *tableForm[] = {(char*)"1I",(char*)"1I",(char*)"1D",(char*)"1D",(char*)"1D",(char*)"1D",(char*)"1D",(char*)"1D",(char*)"1D",(char*)"1D"} ;
    char *tableUnit[] = {(char*)"None",(char*)"None",(char*)"Mpc",(char*)"Mpc",(char*)"Mpc",(char*)"Mpc",(char*)"EeV",(char*)"EeV",(char*)"EeV",(char*)"EeV"} ;
    fits_create_tbl(_fpFitsStream, BINARY_TBL, 0, _fNbFields, tableType, tableForm, tableUnit, extname, &_fStatus) ;
  } else throw TCrpErr("Output format determination failed") ;
  
  if (_fStatus) throw TFitsErr(_fStatus) ;
  _fCurrentRow = 0 ;
  
  // Photon secondary stream if necessary
  if (lpXmlCheckPhoton || lpXmlCheckPPPhoton || (lInteractionType == "Photon")) {
    _fShowerStatus = 0 ;
    fits_open_file(&_fpShowerFitsStream, (aOutputFile).c_str(), READWRITE, &_fShowerStatus) ;
    if (_fShowerStatus) throw TFitsErr(_fShowerStatus) ;
    if (lUnivType == "One Dimension") {
      _fShowerNbFields = 7;
      char extnameShower[] = "PHOTON_SPECTRUM_TABLE" ;
      char *tableTypeShower[] = {(char*)"Particle_Type",
				 (char*)"Cascade_Origin",
				 (char*)"Source_Position",
				 (char*)"Initial_Position",
				 (char*)"Energy",
				 (char*)"Source_Energy",
				 // (char*)"Spectrum_Energy_Bins",
				 (char*)"Spectrum"} ;
      char *tableFormShower[] = {(char*)"1I",(char*)"10A",(char*)"1D",(char*)"1D",(char*)"1D",(char*)"1D",(char*)"170D"} ;
      char *tableUnitShower[] = {(char*)"None",(char*)"None",(char*)"Mpc",(char*)"Mpc",(char*)"EeV",
				 (char*)"EeV",(char*)"Nb of events per log. bin"} ;
      fits_create_tbl(_fpShowerFitsStream, BINARY_TBL, 0, _fShowerNbFields, 
		      tableTypeShower, tableFormShower, tableUnitShower, 
		      extnameShower, &_fShowerStatus) ;
    } else {
      _fShowerNbFields = 17 ;
      char extnameShower[] = "PHOTON_SPECTRUM_TABLE" ;
      char *tableTypeShower[] = {(char*)"Particle_Type",
				 (char*)"Cascade_Origin",
				 (char*)"Initial_Type",
				 (char*)"Source_Position_X",
				 (char*)"Source_Position_Y",
				 (char*)"Source_Position_Z",
				 (char*)"Initial_Position_X",
				 (char*)"Initial_Position_Y",
				 (char*)"Initial_Position_Z",
				 (char*)"Position_X",
				 (char*)"Position_Y",
				 (char*)"Position_Z",
				 (char*)"Energy",
				 (char*)"Momentum_Theta",
				 (char*)"Momentum_Phi",
				 (char*)"Source_Energy",
				 // (char*)"Spectrum_Energy_Bins",
				 (char*)"Spectrum"} ;
      char *tableFormShower[] = {(char*)"1I",(char*)"10A",(char*)"1I",(char*)"1D",(char*)"1D",(char*)"1D",
				 (char*)"1D",(char*)"1D",(char*)"1D",
				 (char*)"1D",(char*)"1D",(char*)"1D",(char*)"1D",
				 (char*)"1D",(char*)"1D",(char*)"1D",(char*)"170D"} ;
      char *tableUnitShower[] = {(char*)"None",(char*)"None",(char*)"None",(char*)"Mpc",(char*)"Mpc",
				 (char*)"Mpc",(char*)"Mpc",(char*)"Mpc",
				 (char*)"Mpc",(char*)"Mpc",(char*)"Mpc",(char*)"Mpc",
				 (char*)"EeV",(char*)"Rad",
				 (char*)"Rad",(char*)"EeV",(char*)"Nb of events per log. bin"};
      fits_create_tbl(_fpShowerFitsStream, BINARY_TBL, 0, _fShowerNbFields, 
		      tableTypeShower, tableFormShower, tableUnitShower, 
		      extnameShower, &_fShowerStatus) ;
    }
    if (_fShowerStatus) throw TFitsErr(_fShowerStatus) ;
    

    _fShowerCurrentRow = 0 ;
    double lEBins[170] = { 1.0066347e-11, 1.2672780e-11, 1.5954085e-11, 2.0085003e-11,
			   2.5285520e-11, 3.1832584e-11, 4.0074849e-11, 5.0451246e-11,
			   6.3514355e-11, 7.9959836e-11, 1.0066347e-10, 1.2672780e-10,
			   1.5954085e-10, 2.0085003e-10, 2.5285520e-10, 3.1832584e-10,
			   4.0074849e-10, 5.0451246e-10, 6.3514355e-10, 7.9959836e-10,
			   1.0066347e-09, 1.2672780e-09, 1.5954085e-09, 2.0085003e-09,
			   2.5285520e-09, 3.1832584e-09, 4.0074849e-09, 5.0451246e-09,
			   6.3514355e-09, 7.9959836e-09, 1.0066347e-08, 1.2672780e-08,
			   1.5954085e-08, 2.0085003e-08, 2.5285520e-08, 3.1832584e-08,
			   4.0074849e-08, 5.0451246e-08, 6.3514355e-08, 7.9959836e-08,
			   1.0066347e-07, 1.2672780e-07, 1.5954085e-07, 2.0085003e-07,
			   2.5285520e-07, 3.1832584e-07, 4.0074849e-07, 5.0451246e-07,
			   6.3514355e-07, 7.9959836e-07, 1.0066347e-06, 1.2672780e-06,
			   1.5954085e-06, 2.0085003e-06, 2.5285520e-06, 3.1832584e-06,
			   4.0074849e-06, 5.0451246e-06, 6.3514355e-06, 7.9959836e-06,
			   1.0066347e-05, 1.2672780e-05, 1.5954085e-05, 2.0085003e-05,
			   2.5285520e-05, 3.1832584e-05, 4.0074849e-05, 5.0451246e-05,
			   6.3514355e-05, 7.9959836e-05, 0.00010066347, 0.00012672780,
			   0.00015954085, 0.00020085003, 0.00025285520, 0.00031832584,
			   0.00040074849, 0.00050451246, 0.00063514355, 0.00079959836,
			   0.0010066347,  0.0012672780,  0.0015954085,  0.0020085003,
			   0.0025285520,  0.0031832584,  0.0040074849,  0.0050451246,
			   0.0063514355,  0.0079959836,  0.010066347,   0.012672780,
			   0.015954085,   0.020085003,   0.025285520,   0.031832584,
			   0.040074849,   0.050451246,   0.063514355,   0.079959836,
			   0.10066347,    0.12672780,    0.15954085,    0.20085003,
			   0.25285520,    0.31832584,    0.40074849,    0.50451246,
			   0.63514355,    0.79959836,    1.0066347,     1.2672780,
			   1.5954085,     2.0085003,     2.5285520,     3.1832584,
			   4.0074849,     5.0451246,     6.3514355,     7.9959836,
			   10.066347,     12.672780,     15.954085,     20.085003,
			   25.285520,     31.832584,     40.074849,     50.451246,
			   63.514355,     79.959836,     100.66347,     126.72780,
			   159.54085,     200.85003,     252.85520,     318.32584,
			   400.74849,     504.51246,     635.14355,     799.59836,
			   1006.6347,     1267.2780,     1595.4085,     2008.5003,
			   2528.5520,     3183.2584,     4007.4849,     5045.1246,
			   6351.4355,     7995.9836,     10066.347,     12672.780,
			   15954.085,     20085.003,     25285.520,     31832.584,
			   40074.849,     50451.246,     63514.355,     79959.836,
			   100663.47,     126727.80,     159540.85,     200850.03,
			   252855.20,     318325.84,     400748.49,     504512.46,
			   635143.55,     799598.36} ; // in EeV
    vector<double> lEnergyBins ;
    for (int i=0; i<170; i++) lEnergyBins.push_back(lEBins[i]) ;
    if (lUnivType == "One Dimension") {
      this->Add1DShower("Ebins_EeV", 0, 0, 0, 0, lEnergyBins) ;
    } else this->Add3DShower("Ebins_EeV", TVector3D(), TVector3D(), TVector3D(), TVector3D(), 0, lEnergyBins, 0) ;
    
  } else _fShowerStatus = -1 ; // flag to say that file was not open
  
  // Neutrino secondary stream if necessary
  if (lpXmlCheckNu) {
    _fNuStatus = 0 ;
    fits_open_file(&_fpNuFitsStream, (aOutputFile).c_str(), READWRITE, &_fNuStatus) ;
    if (_fNuStatus) throw TFitsErr(_fNuStatus) ;
    if (lUnivType == "One Dimension") {
      _fNuNbFields = 7;
      char extnameNu[] = "NEUTRINO_TABLE" ;
      char *tableTypeNu[] = {(char*)"Particle_Type",
			     (char*)"Initial_Type",
			     (char*)"Source_Position",
			     (char*)"Initial_Position",
			     (char*)"Initial_Energy",
			     (char*)"Energy",
			     (char*)"Source_Energy"} ;
      char *tableFormNu[] = {(char*)"1I",(char*)"1I",(char*)"1D",(char*)"1D",(char*)"1D",(char*)"1D",(char*)"1D"} ;
      char *tableUnitNu[] = {(char*)"None",(char*)"None",(char*)"Mpc",(char*)"Mpc",(char*)"EeV",(char*)"EeV",(char*)"EeV"} ;
      fits_create_tbl(_fpNuFitsStream, BINARY_TBL, 0, _fNuNbFields, 
		      tableTypeNu, tableFormNu, tableUnitNu, 
		      extnameNu, &_fNuStatus) ;
    } else {
      _fNuNbFields = 15 ;
      char extnameNu[] = "NEUTRINO_TABLE" ;
      char *tableTypeNu[] = {(char*)"Particle_Type",
			     (char*)"Initial_Type",
			     (char*)"Source_Position_X",
			     (char*)"Source_Position_Y",
			     (char*)"Source_Position_Z",
			     (char*)"Initial_Position_X",
			     (char*)"Initial_Position_Y",
			     (char*)"Initial_Position_Z",
			     (char*)"Position_X",
			     (char*)"Position_Y",
			     (char*)"Position_Z",
			     (char*)"Energy",
			     (char*)"Momentum_Theta",
			     (char*)"Momentum_Phi",
			     (char*)"Source_Energy"} ;
      char *tableFormNu[] = {(char*)"1I",(char*)"1I",(char*)"1D",(char*)"1D",(char*)"1D",
			     (char*)"1D",(char*)"1D",(char*)"1D",
			     (char*)"1D",(char*)"1D",(char*)"1D",(char*)"1D",
			     (char*)"1D",(char*)"1D",(char*)"1D"} ;
      char *tableUnitNu[] = {(char*)"None",(char*)"None",(char*)"Mpc",(char*)"Mpc",
			     (char*)"Mpc",(char*)"Mpc",(char*)"Mpc",
			     (char*)"Mpc",(char*)"Mpc",(char*)"Mpc",(char*)"Mpc",
			     (char*)"EeV",(char*)"Rad",
			     (char*)"Rad",(char*)"EeV"};
      fits_create_tbl(_fpNuFitsStream, BINARY_TBL, 0, _fNuNbFields, 
		      tableTypeNu, tableFormNu, tableUnitNu, 
		      extnameNu, &_fNuStatus) ;
    }
    if (_fNuStatus) throw TFitsErr(_fNuStatus) ;
    _fNuCurrentRow = 0 ;
  } else _fNuStatus = -1 ; // flag to say that file was not open
  
}

void TFitsOutput::Add1DEvent(int aType, 
			     double aInitPosition, 
			     double aInitRedshift,
			     double aInitEnergy, 
			     double aTime, 
			     double aEnergy, int aInitType) {

  int lType = aType ;
  int linitType = aInitType; 
  double lInitPositionMpc = aInitPosition * invMpc ;
  double lInitEnergyEeV = aInitEnergy * invEeV ;
  double lTimeMpc = aTime * c_light * invMpc ;
  double lEnergyEeV = aEnergy * invEeV ;
  fits_insert_rows(_fpFitsStream, _fCurrentRow, 1, &_fStatus) ;
  _fCurrentRow += 1 ;
  fits_write_col(_fpFitsStream, TINT, 1, _fCurrentRow, 1, 1, &lType, &_fStatus) ;
  fits_write_col(_fpFitsStream, TINT, 2, _fCurrentRow, 1, 1, &linitType, &_fStatus) ;
  fits_write_col(_fpFitsStream, TDOUBLE, 3, _fCurrentRow, 1, 1, &lInitPositionMpc, &_fStatus) ;
  fits_write_col(_fpFitsStream, TDOUBLE, 4, _fCurrentRow, 1, 1, &lInitEnergyEeV, &_fStatus) ;
  fits_write_col(_fpFitsStream, TDOUBLE, 5, _fCurrentRow, 1, 1, &lTimeMpc, &_fStatus) ;
  fits_write_col(_fpFitsStream, TDOUBLE, 6, _fCurrentRow, 1, 1, &lEnergyEeV, &_fStatus) ;
  if (_fStatus) throw TFitsErr(_fStatus) ;

}

void TFitsOutput::Add3DEvent(int aType,
			     TVector3D aInitPosition,
			     TVector3D aInitMomentum,
			     double aTime,
			     TVector3D aPosition,
			     TVector3D aMomentum, int aInitType ) {

  int lType = aType ;
  int linitType = aInitType; 
  double lInitPositionXMpc = aInitPosition.x() * invMpc ;
  double lInitPositionYMpc = aInitPosition.y() * invMpc ;
  double lInitPositionZMpc = aInitPosition.z() * invMpc ;
  double lInitEnergyEeV = aInitMomentum.mag() * c_light * invEeV ;
  double lInitMomentumTheta = aInitMomentum.theta() ;
  double lInitMomentumPhi = aInitMomentum.phi() ;
  double lTimeMpc = aTime * c_light * invMpc ;
  double lPositionXMpc = aPosition.x() * invMpc ;
  double lPositionYMpc = aPosition.y() * invMpc ;
  double lPositionZMpc = aPosition.z() * invMpc ;
  double lEnergyEeV = aMomentum.mag() * c_light * invEeV ;
  double lMomentumTheta = aMomentum.theta() ;
  double lMomentumPhi = aMomentum.phi() ;

  fits_insert_rows(_fpFitsStream, _fCurrentRow, 1, &_fStatus) ;
  _fCurrentRow += 1 ;
  fits_write_col(_fpFitsStream, TINT, 1, _fCurrentRow, 1, 1, &lType, &_fStatus) ;
  fits_write_col(_fpFitsStream, TINT, 2, _fCurrentRow, 1, 1, &linitType, &_fStatus) ;
  fits_write_col(_fpFitsStream, TDOUBLE, 3, _fCurrentRow, 1, 1, &lInitPositionXMpc, &_fStatus) ;
  fits_write_col(_fpFitsStream, TDOUBLE, 4, _fCurrentRow, 1, 1, &lInitPositionYMpc, &_fStatus) ;
  fits_write_col(_fpFitsStream, TDOUBLE, 5, _fCurrentRow, 1, 1, &lInitPositionZMpc, &_fStatus) ;
  fits_write_col(_fpFitsStream, TDOUBLE, 6, _fCurrentRow, 1, 1, &lInitEnergyEeV, &_fStatus) ;
  fits_write_col(_fpFitsStream, TDOUBLE, 7, _fCurrentRow, 1, 1, &lInitMomentumTheta, &_fStatus) ;
  fits_write_col(_fpFitsStream, TDOUBLE, 8, _fCurrentRow, 1, 1, &lInitMomentumPhi, &_fStatus) ;
  fits_write_col(_fpFitsStream, TDOUBLE, 9, _fCurrentRow, 1, 1, &lTimeMpc, &_fStatus) ;
  fits_write_col(_fpFitsStream, TDOUBLE, 10, _fCurrentRow, 1, 1, &lPositionXMpc, &_fStatus) ;
  fits_write_col(_fpFitsStream, TDOUBLE, 11, _fCurrentRow, 1, 1, &lPositionYMpc, &_fStatus) ;
  fits_write_col(_fpFitsStream, TDOUBLE, 12, _fCurrentRow, 1, 1, &lPositionZMpc, &_fStatus) ;
  fits_write_col(_fpFitsStream, TDOUBLE, 13, _fCurrentRow, 1, 1, &lEnergyEeV, &_fStatus) ;
  fits_write_col(_fpFitsStream, TDOUBLE, 14, _fCurrentRow, 1, 1, &lMomentumTheta, &_fStatus) ;
  fits_write_col(_fpFitsStream, TDOUBLE, 15, _fCurrentRow, 1, 1, &lMomentumPhi, &_fStatus) ;
  if (_fStatus) throw TFitsErr(_fStatus) ;

}


void TFitsOutput::Add1DTraj(int aType,
			    double aTime,
			    double aPosition,
			    double aEnergy, int aInitType) {
  
  int lType = aType ;
  int linitType=aInitType; 
  double lTimeMpc = aTime * c_light * invMpc ;
  double lPositionMpc = aPosition * invMpc ;
  double lEnergyEeV = aEnergy * invEeV ;

  fits_insert_rows(_fpFitsStream, _fCurrentRow, 1, &_fStatus) ;
  _fCurrentRow += 1 ;
  fits_write_col(_fpFitsStream, TINT, 1, _fCurrentRow, 1, 1, &lType, &_fStatus) ;
  fits_write_col(_fpFitsStream, TINT, 2, _fCurrentRow, 1, 1, &linitType, &_fStatus) ;
  fits_write_col(_fpFitsStream, TDOUBLE, 3, _fCurrentRow, 1, 1, &lTimeMpc, &_fStatus) ;
  fits_write_col(_fpFitsStream, TDOUBLE, 4, _fCurrentRow, 1, 1, &lPositionMpc, &_fStatus) ;
  fits_write_col(_fpFitsStream, TDOUBLE, 5, _fCurrentRow, 1, 1, &lEnergyEeV, &_fStatus) ;
  if (_fStatus) throw TFitsErr(_fStatus) ;

}

void TFitsOutput::Add3DTraj(int aType,
			    double aTime,
			    TVector3D aPosition,
			    TVector3D aMomentum,
			    double aEnergy, int aInitType) {
 
  int lType = aType ;
  int linitType=aInitType; 
  double lTimeMpc = aTime * c_light * invMpc ;
  double lPositionXMpc = aPosition.x() * invMpc ;
  double lPositionYMpc = aPosition.y() * invMpc ;
  double lPositionZMpc = aPosition.z() * invMpc ;
  double lMomentumXEeV = aMomentum.x() * invEeV * c_light;
  double lMomentumYEeV = aMomentum.y() * invEeV * c_light;
  double lMomentumZEeV = aMomentum.z() * invEeV * c_light;
  double lEnergyEeV = aEnergy * invEeV ;

  fits_insert_rows(_fpFitsStream, _fCurrentRow, 1, &_fStatus) ;
  _fCurrentRow += 1 ;
  fits_write_col(_fpFitsStream, TINT, 1, _fCurrentRow, 1, 1, &lType, &_fStatus) ;
  fits_write_col(_fpFitsStream, TINT, 2, _fCurrentRow, 1, 1, &linitType, &_fStatus) ;
  fits_write_col(_fpFitsStream, TDOUBLE, 3, _fCurrentRow, 1, 1, &lTimeMpc, &_fStatus) ;
  fits_write_col(_fpFitsStream, TDOUBLE, 4, _fCurrentRow, 1, 1, &lPositionXMpc, &_fStatus) ;
  fits_write_col(_fpFitsStream, TDOUBLE, 5, _fCurrentRow, 1, 1, &lPositionYMpc, &_fStatus) ;
  fits_write_col(_fpFitsStream, TDOUBLE, 6, _fCurrentRow, 1, 1, &lPositionZMpc, &_fStatus) ;
  fits_write_col(_fpFitsStream, TDOUBLE, 7, _fCurrentRow, 1, 1, &lMomentumXEeV, &_fStatus) ;
  fits_write_col(_fpFitsStream, TDOUBLE, 8, _fCurrentRow, 1, 1, &lMomentumYEeV, &_fStatus) ;
  fits_write_col(_fpFitsStream, TDOUBLE, 9, _fCurrentRow, 1, 1, &lMomentumZEeV, &_fStatus) ;
  fits_write_col(_fpFitsStream, TDOUBLE, 10, _fCurrentRow, 1, 1, &lEnergyEeV, &_fStatus) ;
  if (_fStatus) throw TFitsErr(_fStatus) ;

}

void TFitsOutput::Add3DTraj(int aType,
			    double aTime,
			    TVector3D aPosition,
			    TVector3D aMomentum,
			    double aEnergy,
			    double aMagField, // not implemented in conf...
			    double aDensity, int aInitType) {
 
  int lType = aType ;
  int linitType=aInitType;
  double lTimeMpc = aTime * c_light * invMpc ;
  double lPositionXMpc = aPosition.x() * invMpc ;
  double lPositionYMpc = aPosition.y() * invMpc ;
  double lPositionZMpc = aPosition.z() * invMpc ;
  double lEnergyEeV = aEnergy * invEeV ;
  double lMagFieldMuG = aMagField / muG ;

  fits_insert_rows(_fpFitsStream, _fCurrentRow, 1, &_fStatus) ;
  _fCurrentRow += 1 ;
  fits_write_col(_fpFitsStream, TINT, 1, _fCurrentRow, 1, 1, &lType, &_fStatus) ;
  fits_write_col(_fpFitsStream, TINT, 2, _fCurrentRow, 1, 1, &linitType, &_fStatus) ;
  fits_write_col(_fpFitsStream, TDOUBLE, 3, _fCurrentRow, 1, 1, &lTimeMpc, &_fStatus) ;
  fits_write_col(_fpFitsStream, TDOUBLE, 4, _fCurrentRow, 1, 1, &lPositionXMpc, &_fStatus) ;
  fits_write_col(_fpFitsStream, TDOUBLE, 5, _fCurrentRow, 1, 1, &lPositionYMpc, &_fStatus) ;
  fits_write_col(_fpFitsStream, TDOUBLE, 6, _fCurrentRow, 1, 1, &lPositionZMpc, &_fStatus) ;
  fits_write_col(_fpFitsStream, TDOUBLE, 7, _fCurrentRow, 1, 1, &lEnergyEeV, &_fStatus) ;
  fits_write_col(_fpFitsStream, TDOUBLE, 7, _fCurrentRow, 1, 1, &lMagFieldMuG, &_fStatus) ;
  fits_write_col(_fpFitsStream, TDOUBLE, 7, _fCurrentRow, 1, 1, &aDensity, &_fStatus) ;
  if (_fStatus) throw TFitsErr(_fStatus) ;

}

void TFitsOutput::Add1DShower(string aOrigin,
			      double aSourcePosition,
			      double aInitPosition, 
			      double aEnergy, 
			      double aSourceEnergy,
			      // vector<double> aEnergyBins,
			      vector<double> aSpectrum) {

  int lType = 22 ;
  const char* lOrigin = aOrigin.c_str() ;
  double lSourcePosMpc = aSourcePosition * invMpc ;
  double lInitPosMpc = aInitPosition * invMpc ;
  double lEnergyEeV = aEnergy * invEeV ;
  double lSourceEnergyEeV = aSourceEnergy * invEeV ;
  //  double* lEBins = new double[aEnergyBins.size()] ;
  //for (unsigned int i=0; i<aEnergyBins.size(); i++) lEBins[i] = aEnergyBins[i] * invEeV ;
  double* lSpec = new double[aSpectrum.size()] ;
  for (unsigned int i=0; i<aSpectrum.size(); i++) lSpec[i] = aSpectrum[i] ;
  fits_insert_rows(_fpShowerFitsStream, _fShowerCurrentRow, 1, &_fShowerStatus) ;
  _fShowerCurrentRow += 1 ;
  fits_write_col(_fpShowerFitsStream, TINT, 1, _fShowerCurrentRow, 1, 1, &lType, &_fShowerStatus) ;
  fits_write_col(_fpShowerFitsStream, TSTRING, 2, _fShowerCurrentRow, 1, 1, &lOrigin, &_fShowerStatus) ;
  fits_write_col(_fpShowerFitsStream, TDOUBLE, 3, _fShowerCurrentRow, 1, 1, &lSourcePosMpc, &_fShowerStatus) ;
  fits_write_col(_fpShowerFitsStream, TDOUBLE, 4, _fShowerCurrentRow, 1, 1, &lInitPosMpc, &_fShowerStatus) ;
  fits_write_col(_fpShowerFitsStream, TDOUBLE, 5, _fShowerCurrentRow, 1, 1, &lEnergyEeV, &_fShowerStatus) ;
  fits_write_col(_fpShowerFitsStream, TDOUBLE, 6, _fShowerCurrentRow, 1, 1, &lSourceEnergyEeV, &_fShowerStatus) ;
  fits_write_col(_fpShowerFitsStream, TDOUBLE, 7, _fShowerCurrentRow, 1, 170, lSpec, &_fShowerStatus) ;
  delete[] lSpec ;
  if (_fShowerStatus) throw TFitsErr(_fShowerStatus) ;
  
}

void TFitsOutput::Add3DShower(string aOrigin,
			      TVector3D aSourcePosition,
			      TVector3D aInitPosition, 
			      TVector3D aPosition,
			      TVector3D aMomentum,
			      double aSourceEnergy,
			      //vector<double> aEnergyBins,
			      vector<double> aSpectrum, int aInitType) {

  int lType = 22 ;
  int linitType=aInitType;
  const char* lOrigin = aOrigin.c_str() ;
  double lSourcePosXMpc = aSourcePosition.x() * invMpc ;
  double lSourcePosYMpc = aSourcePosition.y() * invMpc ;
  double lSourcePosZMpc = aSourcePosition.z() * invMpc ;
  double lInitPosXMpc = aInitPosition.x() * invMpc ;
  double lInitPosYMpc = aInitPosition.y() * invMpc ;
  double lInitPosZMpc = aInitPosition.z() * invMpc ;
  double lPosXMpc = aPosition.x() * invMpc ;
  double lPosYMpc = aPosition.y() * invMpc ;
  double lPosZMpc = aPosition.z() * invMpc ;
  double lEnergyEeV = aMomentum.mag() * c_light * invEeV ;
  double lMomentumTheta = aMomentum.theta() ;
  double lMomentumPhi = aMomentum.phi() ;
  double lSourceEnergyEeV = aSourceEnergy * invEeV ;
  //double* lEBins = new double[aEnergyBins.size()] ;
  //for (unsigned int i=0; i<aEnergyBins.size(); i++) lEBins[i] = aEnergyBins[i] * invEeV ;
  double* lSpec = new double[aSpectrum.size()] ;
  for (unsigned int i=0; i<aSpectrum.size(); i++) lSpec[i] = aSpectrum[i] ;

  fits_insert_rows(_fpShowerFitsStream, _fShowerCurrentRow, 1, &_fShowerStatus) ;
  _fShowerCurrentRow += 1 ;
  fits_write_col(_fpShowerFitsStream, TINT, 1, _fShowerCurrentRow, 1, 1, &lType, &_fShowerStatus) ;
  fits_write_col(_fpShowerFitsStream, TSTRING, 2, _fShowerCurrentRow, 1, 1, &lOrigin, &_fShowerStatus) ;
  fits_write_col(_fpShowerFitsStream, TINT, 3, _fShowerCurrentRow, 1, 1, &linitType, &_fShowerStatus) ;
  fits_write_col(_fpShowerFitsStream, TDOUBLE, 4, _fShowerCurrentRow, 1, 1, &lSourcePosXMpc, &_fShowerStatus) ;
  fits_write_col(_fpShowerFitsStream, TDOUBLE, 5, _fShowerCurrentRow, 1, 1, &lSourcePosYMpc, &_fShowerStatus) ;
  fits_write_col(_fpShowerFitsStream, TDOUBLE, 6, _fShowerCurrentRow, 1, 1, &lSourcePosZMpc, &_fShowerStatus) ;
  fits_write_col(_fpShowerFitsStream, TDOUBLE, 7, _fShowerCurrentRow, 1, 1, &lInitPosXMpc, &_fShowerStatus) ;
  fits_write_col(_fpShowerFitsStream, TDOUBLE, 8, _fShowerCurrentRow, 1, 1, &lInitPosYMpc, &_fShowerStatus) ;
  fits_write_col(_fpShowerFitsStream, TDOUBLE, 9, _fShowerCurrentRow, 1, 1, &lInitPosZMpc, &_fShowerStatus) ;
  fits_write_col(_fpShowerFitsStream, TDOUBLE, 10, _fShowerCurrentRow, 1, 1, &lPosXMpc, &_fShowerStatus) ;
  fits_write_col(_fpShowerFitsStream, TDOUBLE, 11, _fShowerCurrentRow, 1, 1, &lPosYMpc, &_fShowerStatus) ;
  fits_write_col(_fpShowerFitsStream, TDOUBLE, 12, _fShowerCurrentRow, 1, 1, &lPosZMpc, &_fShowerStatus) ;
  fits_write_col(_fpShowerFitsStream, TDOUBLE, 13, _fShowerCurrentRow, 1, 1, &lEnergyEeV, &_fShowerStatus) ;
  fits_write_col(_fpShowerFitsStream, TDOUBLE, 14, _fShowerCurrentRow, 1, 1, &lMomentumTheta, &_fShowerStatus) ;
  fits_write_col(_fpShowerFitsStream, TDOUBLE, 15, _fShowerCurrentRow, 1, 1, &lMomentumPhi, &_fShowerStatus) ;
  fits_write_col(_fpShowerFitsStream, TDOUBLE, 16, _fShowerCurrentRow, 1, 1, &lSourceEnergyEeV, &_fShowerStatus) ;
  //fits_write_col(_fpShowerFitsStream, TDOUBLE, 16, _fShowerCurrentRow, 1, 170, lEBins, &_fShowerStatus) ;
  fits_write_col(_fpShowerFitsStream, TDOUBLE, 17, _fShowerCurrentRow, 1, 170, lSpec, &_fShowerStatus) ;
  //  delete[] lEBins;
  delete[] lSpec ;
  if (_fShowerStatus) throw TFitsErr(_fShowerStatus) ;
}


void TFitsOutput::Add1DNeutrino(int aType,
				double aSourcePosition,
				double aInitPosition, 
				double aInitEnergy,
				double aEnergy,
				double aSourceEnergy, int aInitType) {

  int lType = aType ;
  int linitType= aInitType;
  double lSourcePosMpc = aSourcePosition * invMpc ;
  double lInitPosMpc = aInitPosition * invMpc ;
  double lInitEnergyEeV = aInitEnergy * invEeV ;
  double lEnergyEeV = aEnergy * invEeV ;
  double lSourceEnergyEeV = aSourceEnergy * invEeV ;
  fits_insert_rows(_fpNuFitsStream, _fNuCurrentRow, 1, &_fNuStatus) ;
  _fNuCurrentRow += 1 ;
  fits_write_col(_fpNuFitsStream, TINT, 1, _fNuCurrentRow, 1, 1, &lType, &_fNuStatus) ;
  fits_write_col(_fpNuFitsStream, TINT, 2, _fNuCurrentRow, 1, 1, &linitType, &_fNuStatus) ;
  fits_write_col(_fpNuFitsStream, TDOUBLE, 3, _fNuCurrentRow, 1, 1, &lSourcePosMpc, &_fNuStatus) ;
  fits_write_col(_fpNuFitsStream, TDOUBLE, 4, _fNuCurrentRow, 1, 1, &lInitPosMpc, &_fNuStatus) ;
  fits_write_col(_fpNuFitsStream, TDOUBLE, 5, _fNuCurrentRow, 1, 1, &lInitEnergyEeV, &_fNuStatus) ;
  fits_write_col(_fpNuFitsStream, TDOUBLE, 6, _fNuCurrentRow, 1, 1, &lEnergyEeV, &_fNuStatus) ;
  fits_write_col(_fpNuFitsStream, TDOUBLE, 7, _fNuCurrentRow, 1, 1, &lSourceEnergyEeV, &_fNuStatus) ;
  if (_fNuStatus) throw TFitsErr(_fNuStatus) ;
}


void TFitsOutput::Add3DNeutrino(int aType,
				TVector3D aSourcePosition,
				TVector3D aInitPosition, 
				TVector3D aPosition,
				TVector3D aMomentum,
				double aSourceEnergy, int aInitType) {

  int lType = aType ;
  int linitType= aInitType; 
  double lSourcePosXMpc = aSourcePosition.x() * invMpc ;
  double lSourcePosYMpc = aSourcePosition.y() * invMpc ;
  double lSourcePosZMpc = aSourcePosition.z() * invMpc ;
  double lInitPosXMpc = aInitPosition.x() * invMpc ;
  double lInitPosYMpc = aInitPosition.y() * invMpc ;
  double lInitPosZMpc = aInitPosition.z() * invMpc ;
  double lPosXMpc = aPosition.x() * invMpc ;
  double lPosYMpc = aPosition.y() * invMpc ;
  double lPosZMpc = aPosition.z() * invMpc ;
  double lEnergyEeV = aMomentum.mag() * c_light * invEeV ;
  double lMomentumTheta = aMomentum.theta() ;
  double lMomentumPhi = aMomentum.phi() ;
  double lSourceEnergyEeV = aSourceEnergy * invEeV ;

  fits_insert_rows(_fpNuFitsStream, _fNuCurrentRow, 1, &_fNuStatus) ;
  _fNuCurrentRow += 1 ;
  fits_write_col(_fpNuFitsStream, TINT, 1, _fNuCurrentRow, 1, 1, &lType, &_fNuStatus) ;
  fits_write_col(_fpNuFitsStream, TINT, 2, _fNuCurrentRow, 1, 1, &linitType, &_fNuStatus) ;
  fits_write_col(_fpNuFitsStream, TDOUBLE, 3, _fNuCurrentRow, 1, 1, &lSourcePosXMpc, &_fNuStatus) ;
  fits_write_col(_fpNuFitsStream, TDOUBLE, 4, _fNuCurrentRow, 1, 1, &lSourcePosYMpc, &_fNuStatus) ;
  fits_write_col(_fpNuFitsStream, TDOUBLE, 5, _fNuCurrentRow, 1, 1, &lSourcePosZMpc, &_fNuStatus) ;
  fits_write_col(_fpNuFitsStream, TDOUBLE, 6, _fNuCurrentRow, 1, 1, &lInitPosXMpc, &_fNuStatus) ;
  fits_write_col(_fpNuFitsStream, TDOUBLE, 7, _fNuCurrentRow, 1, 1, &lInitPosYMpc, &_fNuStatus) ;
  fits_write_col(_fpNuFitsStream, TDOUBLE, 8, _fNuCurrentRow, 1, 1, &lInitPosZMpc, &_fNuStatus) ;
  fits_write_col(_fpNuFitsStream, TDOUBLE, 9, _fNuCurrentRow, 1, 1, &lPosXMpc, &_fNuStatus) ;
  fits_write_col(_fpNuFitsStream, TDOUBLE, 10, _fNuCurrentRow, 1, 1, &lPosYMpc, &_fNuStatus) ;
  fits_write_col(_fpNuFitsStream, TDOUBLE, 11, _fNuCurrentRow, 1, 1, &lPosZMpc, &_fNuStatus) ;
  fits_write_col(_fpNuFitsStream, TDOUBLE, 12, _fNuCurrentRow, 1, 1, &lEnergyEeV, &_fNuStatus) ;
  fits_write_col(_fpNuFitsStream, TDOUBLE, 13, _fNuCurrentRow, 1, 1, &lMomentumTheta, &_fNuStatus) ;
  fits_write_col(_fpNuFitsStream, TDOUBLE, 14, _fNuCurrentRow, 1, 1, &lMomentumPhi, &_fNuStatus) ;
  fits_write_col(_fpNuFitsStream, TDOUBLE, 15, _fNuCurrentRow, 1, 1, &lSourceEnergyEeV, &_fNuStatus) ;
  if (_fNuStatus) throw TFitsErr(_fNuStatus) ;
}


TFitsOutput::~TFitsOutput() {
  fits_close_file(_fpFitsStream,&_fStatus) ;
  if (_fStatus) throw TFitsErr(_fStatus) ;
  if (_fShowerStatus != -1 ) {
    fits_close_file(_fpShowerFitsStream,&_fShowerStatus) ;
    if (_fShowerStatus) throw TFitsErr(_fShowerStatus) ;
  }
  if (_fNuStatus != -1 ) {
    fits_close_file(_fpNuFitsStream,&_fNuStatus) ;
    if (_fNuStatus) throw TFitsErr(_fNuStatus) ;
  }
}
