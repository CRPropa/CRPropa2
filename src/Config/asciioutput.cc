

/**
   @file   asciioutput.cc
   @author Eric Armengaud, armengau@in2p3.fr
   @brief  Implementation of the TAsciiOutput class. See the .h file
*/

#include "asciioutput.h"

TAsciiOutput::TAsciiOutput(const char* aFileName, string aOutputFile, bool aForceOverwrite ) : TXmlParam(aFileName){
	SetType(OUTPUT_ASCII);
  _fFileName = aOutputFile ;

  string lUnivType=(XmlExtract().GetElement("Environment"))->Attribute("type");
  string lRecordType = (XmlExtract().GetElement("Output"))->Attribute("type");
  TiXmlElement* lpXmlInter = XmlExtract().GetElement("Interactions") ;
  TiXmlElement* lpXmlPhoton = lpXmlInter->FirstChildElement("SecondaryPhotons") ;
  TiXmlElement* lpXmlPPPhoton = lpXmlInter->FirstChildElement("SecondaryPairProdPhotons") ;
  TiXmlElement* lpXmlNu = lpXmlInter->FirstChildElement("SecondaryNeutrinos") ;
  string lInteractionType = lpXmlInter->Attribute("type") ;

  struct stat lFilestat ; // man 2 stat :-)
  if ( stat(aOutputFile.c_str(), &lFilestat) == 0 && ! aForceOverwrite )
    throw TCrpErr("Output file already exists. If you want to overwrite, use option=\"force\".") ;

  _fDataStream.open(aOutputFile.c_str(), ios::out ) ;
  if (_fDataStream.is_open()) {
    // "Header" filling
    _fDataStream << "#CRPropa - Output data file" << endl ;
    if (lUnivType == "One Dimension" && lRecordType == "Events") {
      _fDataStream << "#Format - Particle_Type Initial_Particle_Type Initial_Position(Mpc) Initial_Energy(EeV) Time(Mpc) Energy(EeV)" << endl ;
    } else if (lUnivType == "One Dimension" && lRecordType == "Full Trajectories") {
      _fDataStream << "#Format - Particle_Type Initial_Particle_Type Time(Mpc) Position(Mpc) Energy(EeV)" << endl ;
    } else if (lUnivType != "One Dimension" && lRecordType == "Events") {
      _fDataStream << "#Format - Particle_Type Initial_Particle_Type Initial_Position[X,Y,Z](Mpc) Initial_Momentum[E,theta,phi](EeV) Time(Mpc) Position[X,Y,Z](Mpc) Momentum[E,theta,phi](EeV)" << endl ;
    } else if (lUnivType != "One Dimension" && lRecordType == "Full Trajectories") {
      _fDataStream << "#Format - Particle_Type Initial_Particle_Type Time(Mpc) Position[X,Y,Z](Mpc) Momentum[X,Y,Z](Mpc) Energy(EeV)" << endl ;
    } else throw TCrpErr("Output format determination failed") ;

  } else {
    throw TCrpErr("Ascii output file cannot be opened");
  }
 
  if (!_fDataStream.good()) throw TCrpErr("Error with ascii output") ;

  // Photon secondary stream if necessary
  if (lpXmlPhoton || lpXmlPPPhoton || (lInteractionType == "Photon")) {
    if ( stat((aOutputFile+".photons").c_str(), &lFilestat) == 0 && ! aForceOverwrite )
      throw TCrpErr("Output temporary photon file already exists.") ;
    _fShowerStream.open((aOutputFile+".photons").c_str(), ios::out) ;
    if (_fShowerStream.is_open()) {
      _fShowerStream << "#Output secondary photons" << endl ;
      if (lUnivType == "One Dimension") {
	_fShowerStream << "#Format - Particle_Type Cascade_Origin Init_Type Source_Position(Mpc) Initial_Position(Mpc) Energy(EeV) Source_Energy(EeV) - Spectrum(Nb of events per log. bin)" << endl ;
      } else {
	_fShowerStream << "#Format - Particle_Type Cascade_Origin  Init_Type Source_Position[X,Y,Z](Mpc) Initial_Position[X,Y,Z](Mpc) Position[X,Y,Z](Mpc) Momentum[E,theta,phi](EeV) Source_Energy(EeV) - Spectrum(Nb of events per log. bin)" << endl ;
      }
    } else throw TCrpErr("Ascii shower output file cannot be opened");

    if (!_fShowerStream.good()) throw TCrpErr("Error with ascii output") ;
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
    } else {
      this->Add3DShower("Ebins_EeV", TVector3D(), TVector3D(), TVector3D(), TVector3D(), 0, lEnergyBins, 0) ;
    }
    
  }
  
  // Neutrino secondary stream if necessary
  if (lpXmlNu) {
    if ( stat((aOutputFile+".neutrinos").c_str(), &lFilestat) == 0 && ! aForceOverwrite )
      throw TCrpErr("Output temporary neutrino file already exists.") ;
    _fNeutrinoStream.open((aOutputFile+".neutrinos").c_str(), ios::out) ;
    if (_fNeutrinoStream.is_open()) {
      _fNeutrinoStream << "#Output secondary neutrinos" << endl ;
      if (lUnivType == "One Dimension") {
	_fNeutrinoStream << "#Format - Particle_Type Generating_Particle Source_Position(Mpc) Initial_Position(Mpc) Initial_Energy(EeV) Energy(EeV) Source_Energy(EeV)" << endl ;
      } else {
	_fNeutrinoStream << "#Format - Particle_Type Generating_Particle Source_Position[X,Y,Z](Mpc) Initial_Position[X,Y,Z](Mpc) Position[X,Y,Z](Mpc) Momentum[E,theta,phi](EeV) Source_Energy(EeV)" << endl ;
      }
    } else throw TCrpErr("Ascii neutrino output file cannot be opened");
    
    if (!_fNeutrinoStream.good()) throw TCrpErr("Error with ascii output") ;
  }
  
}

void TAsciiOutput::Add1DEvent(int aType,
			      double aInitPosition, 
			      double aInitRedshift,
			      double aInitEnergy, 
			      double aTime, 
			      double aEnergy, int aInitType)  {
  
  _fDataStream << aType << " " << aInitType << " "
	       << aInitPosition * invMpc << " "
	       << aInitEnergy * invEeV << " "
	       << aTime * c_light * invMpc << " "
	       << aEnergy * invEeV << endl ;
  if (!_fDataStream.good()) throw TCrpErr("Error with ascii output") ;

}

void TAsciiOutput::Add3DEvent(int aType,
			      TVector3D aInitPosition,
			      TVector3D aInitMomentum,
			      double aTime,
			      TVector3D aPosition,
			      TVector3D aMomentum,  int aInitType) {

  _fDataStream << aType << " "<< aInitType << " "
	       << aInitPosition.x() * invMpc << " "
	       << aInitPosition.y() * invMpc << " "
	       << aInitPosition.z() * invMpc << " "
	       << aInitMomentum.mag() * c_light * invEeV << " "
	       << aInitMomentum.theta() << " "
	       << aInitMomentum.phi() << " "
	       << aTime * c_light * invMpc << " "
	       << aPosition.x() * invMpc << " "
	       << aPosition.y() * invMpc << " "
	       << aPosition.z() * invMpc << " "
	       << aMomentum.mag() * c_light * invEeV << " "
	       << aMomentum.theta() << " "
	       << aMomentum.phi() << endl ;
  if (!_fDataStream.good()) throw TCrpErr("Error with ascii output") ;

}


void TAsciiOutput::Add1DTraj(int aType,
			     double aTime,
			     double aPosition,
			     double aEnergy,  int aInitType) {
  
  _fDataStream << aType << " " << aInitType << " " 
	       << aTime * c_light * invMpc << " " // rq: time != pos if secondary parts.(propa)
	       << aPosition * invMpc << " "
	       << aEnergy * invEeV << endl ;
  if (!_fDataStream.good()) throw TCrpErr("Error with ascii output") ;

}

void TAsciiOutput::Add3DTraj(int aType,
			     double aTime,
			     TVector3D aPosition,
			     TVector3D aMomentum,
			     double aEnergy, int aInitType) {

  _fDataStream << aType << " " << aInitType << " " 
	       << aTime * c_light * invMpc << " "
	       << aPosition.x() * invMpc << " "
	       << aPosition.y() * invMpc << " "
	       << aPosition.z() * invMpc << " "
	       << aMomentum.x() * invEeV * c_light << " "
		   << aMomentum.y() * invEeV * c_light << " "
		   << aMomentum.z() * invEeV * c_light << " "
	       << aEnergy * invEeV << endl ;
  if (!_fDataStream.good()) throw TCrpErr("Error with ascii output") ;

}

void TAsciiOutput::Add3DTraj(int aType,
			     double aTime,
			     TVector3D aPosition,
			     TVector3D aMomentum,
			     double aEnergy,
			     double aMagField, // not implemented in conf. files...
			     double aDensity, int aInitType) {

  _fDataStream << aType << " " <<  aInitType << " "
	       << aTime * invMpc << " "
	       << aPosition.x() * invMpc << " "
	       << aPosition.y() * invMpc << " "
	       << aPosition.z() * invMpc << " "
	       << aMomentum.x() * invEeV << " "
		   << aMomentum.y() * invEeV << " "
		   << aMomentum.z() * invEeV << " "
	       << aEnergy * invEeV << " "
	       << aMagField / muG << " "
	       << aDensity << endl ;
  if (!_fDataStream.good()) throw TCrpErr("Error with ascii output") ;

}

void TAsciiOutput::Add1DShower(string aOrigin,
			       double aSourcePosition,
			       double aInitPosition, 
			       double aEnergy, 
			       double aSourceEnergy,
			       // vector<double> aEnergyBins,
			       vector<double> aSpectrum) {

  unsigned int lSpectrumSize = aSpectrum.size();
  //if (aEnergyBins.size() != lSpectrumSize) throw TCrpErr("Error: size of spectrum in Add1DShower") ;

  _fShowerStream << 22 << " "
		 << aOrigin << " "
		 << aSourcePosition * invMpc << " "
		 << aInitPosition * invMpc << " "
		 << aEnergy * invEeV << " "
		 << aSourceEnergy * invEeV << endl ;
  //for (unsigned int i=0; i<lSpectrumSize; i++) _fShowerStream << aEnergyBins[i]* invEeV << " " ; // in EeV.
  //_fShowerStream << endl ;
  for (unsigned int i=0; i<lSpectrumSize; i++) _fShowerStream << aSpectrum[i] << " " ;
  _fShowerStream << endl ;
  if (!_fShowerStream.good()) throw TCrpErr("Error with ascii output") ;

}

void TAsciiOutput::Add3DShower(string aOrigin,
			       TVector3D aSourcePosition,
			       TVector3D aInitPosition, 
			       TVector3D aPosition,
			       TVector3D aMomentum,
			       double aSourceEnergy,
			       // vector<double> aEnergyBins,
			       vector<double> aSpectrum,  int aInitType) {

  unsigned int lSpectrumSize = aSpectrum.size();
  //if (aEnergyBins.size() != lSpectrumSize) throw TCrpErr("pbl size of spectrum in add3dshower") ;

  _fShowerStream << 22 << " "
		 << aOrigin << " " << aInitType << " " 
		 << aSourcePosition.x() * invMpc << " "
		 << aSourcePosition.y() * invMpc << " "
		 << aSourcePosition.z() * invMpc << " "
		 << aInitPosition.x() * invMpc << " "
		 << aInitPosition.y() * invMpc << " "
		 << aInitPosition.z() * invMpc << " "
		 << aPosition.x() * invMpc << " "
		 << aPosition.y() * invMpc << " "
		 << aPosition.z() * invMpc << " "
		 << aMomentum.mag() * c_light * invEeV << " "
		 << aMomentum.theta() << " "
		 << aMomentum.phi() << " "
		 << aSourceEnergy * invEeV << endl ;
  //for (unsigned int i=0; i<lSpectrumSize; i++) _fShowerStream << aEnergyBins[i] * invEeV << " " ; // EeV.
  //_fShowerStream << endl ;
  for (unsigned int i=0; i<lSpectrumSize; i++) _fShowerStream << aSpectrum[i] << " " ; // units not clear!
  _fShowerStream << endl ;
  if (!_fShowerStream.good()) throw TCrpErr("Error with ascii output") ;

}

void TAsciiOutput::Add1DNeutrino(int aType,
				 double aSourcePosition,
				 double aInitPosition,
				 double aInitEnergy, 
				 double aEnergy,
				 double aSourceEnergy, int aInitType) {

  _fNeutrinoStream << aType << " " << aInitType << " "
		   << aSourcePosition * invMpc << " "
		   << aInitPosition * invMpc << " "
		   << aInitEnergy * invEeV << " "
		   << aEnergy * invEeV << " "
		   << aSourceEnergy * invEeV << endl ;
  if (!_fNeutrinoStream.good()) throw TCrpErr("Error with ascii output") ;

}

void TAsciiOutput::Add3DNeutrino(int aType,
				 TVector3D aSourcePosition,
				 TVector3D aInitPosition, 
				 TVector3D aPosition,
				 TVector3D aMomentum,
				 double aSourceEnergy,  int aInitType) {

  _fNeutrinoStream << aType << " " << aInitType << " "
		   << aSourcePosition.x() * invMpc << " "
		   << aSourcePosition.y() * invMpc << " "
		   << aSourcePosition.z() * invMpc << " "
		   << aInitPosition.x() * invMpc << " "
		   << aInitPosition.y() * invMpc << " "
		   << aInitPosition.z() * invMpc << " "
		   << aPosition.x() * invMpc << " "
		   << aPosition.y() * invMpc << " "
		   << aPosition.z() * invMpc << " "
		   << aMomentum.mag() * c_light * invEeV << " "
		   << aMomentum.theta() << " "
		   << aMomentum.phi() << " "
		   << aSourceEnergy * invEeV << endl ;
  if (!_fNeutrinoStream.good()) throw TCrpErr("Error with ascii output") ;

}


TAsciiOutput::~TAsciiOutput() {
  _fDataStream.close() ;
  if (_fShowerStream.is_open()) {
    _fShowerStream.close() ;
    string lCommand1 = "cat "+_fFileName+" "+_fFileName+".photons > "+_fFileName+".tmp" ;
    string lCommand2 = "mv -f "+_fFileName+".tmp "+_fFileName ;
    string lCommand3 = "rm -f "+_fFileName+".photons" ;
    int lCheck1 = system(lCommand1.c_str()) ;
    int lCheck2 = system(lCommand2.c_str()) ;
    int lCheck3 = system(lCommand3.c_str()) ;
    if (lCheck1 == -1 || lCheck2 == -1 || lCheck3 == -1) throw TCrpErr("pbl dans mvt des fichiers ascii d'output..") ;
  }
  if (_fNeutrinoStream.is_open()) {
    _fNeutrinoStream.close() ;
    string lCommand1 = "cat "+_fFileName+" "+_fFileName+".neutrinos > "+_fFileName+".tmp" ;
    string lCommand2 = "mv -f "+_fFileName+".tmp "+_fFileName ;
    string lCommand3 = "rm -f "+_fFileName+".neutrinos" ;
    int lCheck1 = system(lCommand1.c_str()) ;
    int lCheck2 = system(lCommand2.c_str()) ;
    int lCheck3 = system(lCommand3.c_str()) ;
    if (lCheck1 == -1 || lCheck2 == -1 || lCheck3 == -1) throw TCrpErr("pbl dans mvt des fichiers ascii d'output..") ;
  }
}
