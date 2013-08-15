/**
   @file   rootoutput.cc
   @author Eric Armengaud, armengau@in2p3.fr
   @brief  Implementation of the TRootOutput class. See the .h file
*/

#include "rootoutput.h"
//#ifdef HAVE_TROOT_H
using namespace std;
TRootOutput::TRootOutput(const char* aFileName, string aOutputFile, bool aForceOverwrite ) : TXmlParam(aFileName){
  SetType(OUTPUT_ROOT);
  _fFileName = aOutputFile ;
  string lUnivType=(XmlExtract().GetElement("Environment"))->Attribute("type");
  string lRecordType = (XmlExtract().GetElement("Output"))->Attribute("type");
  TiXmlElement* lpXmlInter = XmlExtract().GetElement("Interactions") ;
  TiXmlElement* lpXmlPhoton = lpXmlInter->FirstChildElement("SecondaryPhotons") ;
  TiXmlElement* lpXmlPPPhoton = lpXmlInter->FirstChildElement("SecondaryPairProdPhotons") ;
  TiXmlElement* lpXmlNu = lpXmlInter->FirstChildElement("SecondaryNeutrinos") ;
  string lInteractionType = lpXmlInter->Attribute("type") ;

  struct stat lFilestat ; // man 2 stat
  if ( stat(aOutputFile.c_str(), &lFilestat) == 0 && ! aForceOverwrite )
    throw TCrpErr("Output file already exists. If you want to overwrite, use option=\"force\".") ;
  _fRootFile = new TFile(aOutputFile.c_str(),"RECREATE","CRPropa output data file");
  if (_fRootFile->IsZombie()) throw TCrpErr("Root output file cannot be properly opened");
  if (lUnivType == "One Dimension" && lRecordType == "Events") {
    _fEventNtuple = new TNtuple("events","CRPropa 1D events","Particle_Type:Initial_Type:Initial_Position_Mpc:Inital_Redshift:Initial_Energy_EeV:Time_Mpc:Energy_EeV");
  } else if (lUnivType == "One Dimension" && lRecordType == "Full Trajectories") {
    _fEventNtuple = new TNtuple("traj","CRPropa 1D trajectories","Particle_Type:Initial_Type:Time_Mpc:Position_Mpc:Energy_EeV");
  } else if (lUnivType != "One Dimension" && lRecordType == "Events") {
    _fEventNtuple = new TNtuple("events","CRPropa 3D events","Particle_Type:Initial_Type:Initial_Position_X_Mpc:Initial_Position_Y_Mpc:Initial_Position_Z_Mpc:Initial_Momentum_E_EeV:Initial_Momentum_theta:Initial_Momentum_phi:Time_Mpc:Position_X_Mpc:Position_Y_Mpc:Position_Z_Mpc:Momentum_E_EeV:Momentum_theta:Momentum_phi");
  } else if (lUnivType != "One Dimension" && lRecordType == "Full Trajectories") {
    _fEventNtuple = new TNtuple("traj","CRPropa 3D trajectories","Particle_Type:Initial_Type:Time_Mpc:Position_X_Mpc:Position_Y_Mpc:Position_Z_Mpc:Momentum_X_EeV:Momentum_Y_EeV:Momentum_Z_EeV:Energy_EeV");
  } else throw TCrpErr("Output format determination failed") ;
  
 // Neutrino secondary stream if necessary
  if (lpXmlNu) {
    if (lUnivType == "One Dimension") {
      _fNeutrinoNtuple = new TNtuple("neutrinos","CRPropa 1D events","Particle_Type:Initial_Type:Source_Position_Mpc:Initial_Position_Mpc:Initial_Energy_EeV:Energy_EeV:Source_Energy_EeV");
    } else {
      _fNeutrinoNtuple = new TNtuple("neutrinos","CRPropa 3D events","Particle_Type:Initial_Type:Source_Position_X_Mpc:Source_Position_Y_Mpc:Source_Position_Z_Mpc:Initial_Position_X_Mpc:Initial_Position_Y_Mpc:Initial_Position_Z_Mpc:Position_X_Mpc:Position_Y_Mpc:Position_Z_Mpc:Momentum_E_EeV:Momentum_theta:Momentum_phi:Source_Energy_EeV");
    }
  }

// Photon secondary stream if necessary
  if (lpXmlPhoton || lpXmlPPPhoton || (lInteractionType == "Photon")) {
    if (lUnivType == "One Dimension") {
      _fShowerTree = new TTree("photons","CRPropa 1D electromagnetic cascades");
      _fShowerTree->Branch("Particle_Type",&_fParticle_Type,"Particle_Type/I");
      _fShowerTree->Branch("Initial_Type",&_fInit_Particle_Type,"Init_Particle_Type/I");
      _fShowerTree->Branch("Cascade_Origin_String",&_fCascadeOriginString,"Cascade_Origin_String/C");
      _fShowerTree->Branch("SourcePosition_Mpc",&_fSourcePosition_X_Mpc,"SourcePosition_Mpc/D");
      _fShowerTree->Branch("InitialPosition_Mpc",&_fInitialPosition_X_Mpc,"InitialPosition_Mpc/D");
    } else {
      _fShowerTree = new TTree("photons","CRPropa 3D electromagnetic cascades");
      _fShowerTree->Branch("Particle_Type",&_fParticle_Type,"Particle_Type/I");
      _fShowerTree->Branch("Initial_Type",&_fInit_Particle_Type,"Init_Particle_Type/I");
      _fShowerTree->Branch("Cascade_Origin_String",&_fCascadeOriginString,"Cascade_Origin_String/C");
      _fShowerTree->Branch("SourcePosition_X_Mpc",&_fSourcePosition_X_Mpc,"SourcePosition_X_Mpc/D");
      _fShowerTree->Branch("InitialPosition_X_Mpc",&_fInitialPosition_X_Mpc,"InitialPosition_X_Mpc/D");
      _fShowerTree->Branch("Position_X_Mpc",&_fPosition_X_Mpc,"Position_X_Mpc/D");
      _fShowerTree->Branch("SourcePosition_Y_Mpc",&_fSourcePosition_Y_Mpc,"SourcePosition_Y_Mpc/D");
      _fShowerTree->Branch("InitialPosition_Y_Mpc",&_fInitialPosition_Y_Mpc,"InitialPosition_Y_Mpc/D");
      _fShowerTree->Branch("Position_Y_Mpc",&_fPosition_Y_Mpc,"Position_Y_Mpc/D");
      _fShowerTree->Branch("SourcePosition_Z_Mpc",&_fSourcePosition_Z_Mpc,"SourcePosition_Z_Mpc/D");
      _fShowerTree->Branch("InitialPosition_Z_Mpc",&_fInitialPosition_Z_Mpc,"InitialPosition_Z_Mpc/D");
      _fShowerTree->Branch("Position_Z_Mpc",&_fPosition_Z_Mpc,"Position_Z_Mpc/D");
      _fShowerTree->Branch("Momentum_theta",&_fMomentum_theta,"Momentum_theta/D");
      _fShowerTree->Branch("Momentum_phi",&_fMomentum_phi,"Momentum_phi/D");
    }
    _fShowerTree->Branch("Energy_EeV",&_fEnergy_EeV,"Energy_EeV/D");
    _fShowerTree->Branch("SourceEnergy_EeV",&_fSourceEnergy_EeV,"SourceEnergy_EeV/D");
    _fShowerTree->Branch("SpectrumEnergy_EeV",_fSpectrumEnergy_EeV,"SpectrumEnergy_EeV[170]/F");
    _fShowerTree->Branch("Spectrum",_fSpectrum,"Spectrum[170]/F");
    float lEBins[170] = { 1.0066347e-11, 1.2672780e-11, 1.5954085e-11, 2.0085003e-11,
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
     for (int k=0; k<170; k++) _fSpectrumEnergy_EeV[k]=lEBins[k];
  }

}

void TRootOutput::Add1DEvent(int aType,
			     double aInitPosition,
			     double aInitRedshift,
			     double aInitEnergy, 
			     double aTime, 
			     double aEnergy, int aInitType)  {
  _fEventNtuple->Fill(aType, aInitType,
		     aInitPosition*invMpc,
		     aInitRedshift, 
		     aInitEnergy*invEeV,
		     aTime*c_light*invMpc,
		     aEnergy*invEeV);
}

void TRootOutput::Add3DEvent(int aType,
			     TVector3D aInitPosition,
			     TVector3D aInitMomentum,
			     double aTime,
			     TVector3D aPosition,
			     TVector3D aMomentum, int aInitType) {
  _fEventNtuple->Fill(aType, aInitType,
		     aInitPosition.x() * invMpc,
		     aInitPosition.y() * invMpc,
		     aInitPosition.z() * invMpc,
		     aInitMomentum.mag() * c_light * invEeV,
		     aInitMomentum.theta(),
		     aInitMomentum.phi(),
		     aTime * c_light * invMpc,
		     aPosition.x() * invMpc,
		     aPosition.y() * invMpc,
		     aPosition.z() * invMpc,
		     aMomentum.mag() * c_light * invEeV,
		     aMomentum.theta(),
		     aMomentum.phi());
}

void TRootOutput::Add1DTraj(int aType,
			    double aTime,
			    double aPosition,
			    double aEnergy, int aInitType) {
  _fEventNtuple->Fill(aType, aInitType,
		      aTime*c_light*invMpc,
		      aPosition*invMpc,
		      aEnergy*invEeV);
}

void TRootOutput::Add3DTraj(int aType,
			    double aTime,
			    TVector3D aPosition,
			    TVector3D aMomentum,
			    double aEnergy, int aInitType) {
  _fEventNtuple->Fill(aType, aInitType,
		      aTime*c_light*invMpc,
		      aPosition.x()*invMpc,
		      aPosition.y()*invMpc,
		      aPosition.z()*invMpc,
		      aMomentum.x()*invEeV* c_light,
			  aMomentum.y()*invEeV* c_light,
			  aMomentum.z()*invEeV* c_light ,
		      aEnergy*invEeV);
}

void TRootOutput::Add3DTraj(int aType,
			    double aTime,
			    TVector3D aPosition,
			    TVector3D aMomentum,
			    double aEnergy,
			    double aMagField,
			    double aDensity, int aInitType) {
  _fEventNtuple->Fill(aType,  aInitType,
		      aTime*c_light*invMpc,
		      aPosition.x()*invMpc,
		      aPosition.y()*invMpc,
		      aPosition.z()*invMpc,
		      aMomentum.x()*invEeV,
			  aMomentum.y()*invEeV,
			  aMomentum.z()*invEeV,
		      aEnergy*invEeV,
		      aMagField/muG,
		      aDensity);
}

void TRootOutput::Add1DShower(string aOrigin, 
			      //int aInitType,
			      double aSourcePosition,
			      double aInitPosition, 
			      double aEnergy, 
			      double aSourceEnergy,
			      // vector<double> aEnergyBins,
			      vector<double> aSpectrum) {
  _fParticle_Type = 22 ;
  sprintf(_fCascadeOriginString, aOrigin.c_str());
  //_fInit_Particle_Type = aInitType;
  //if      (aOrigin == "all")       _fCascadeOrigin=0;
  //else if (aOrigin == "gzk_gamma") _fCascadeOrigin=1;
  //else if (aOrigin == "gzk_e+")    _fCascadeOrigin=2;
  //else if (aOrigin == "gzk_e-")    _fCascadeOrigin=3;
  //else if (aOrigin == "pairprod")  _fCascadeOrigin=4;
  //else throw TCrpErr("Unknown code for the EM cascade origin") ;
  _fSourcePosition_X_Mpc = aSourcePosition*invMpc;
  _fInitialPosition_X_Mpc = aInitPosition*invMpc;
  _fEnergy_EeV = aEnergy*invEeV;
  _fSourceEnergy_EeV = aSourceEnergy*invEeV;
  for (int k=0;k<170;k++) (_fSpectrum)[k]=(float)(aSpectrum[k]);
  _fShowerTree->Fill();
}

void TRootOutput::Add3DShower(string aOrigin,
			       TVector3D aSourcePosition,
			       TVector3D aInitPosition, 
			       TVector3D aPosition,
			       TVector3D aMomentum, 
			       double aSourceEnergy,
			       // vector<double> aEnergyBins,
			       vector<double> aSpectrum, int aInitType) {
  _fParticle_Type = 22 ;
  _fInit_Particle_Type = aInitType;
  sprintf(_fCascadeOriginString, aOrigin.c_str());
  //if      (aOrigin == "all")       _fCascadeOrigin=0;
  //else if (aOrigin == "gzk_gamma") _fCascadeOrigin=1;
  //else if (aOrigin == "gzk_e+")    _fCascadeOrigin=2;
  //else if (aOrigin == "gzk_e-")    _fCascadeOrigin=3;
  //else if (aOrigin == "pairprod")  _fCascadeOrigin=4;
  _fSourcePosition_X_Mpc = aSourcePosition.x()*invMpc;
  _fInitialPosition_X_Mpc = aInitPosition.x()*invMpc;
  _fPosition_X_Mpc = aPosition.x()*invMpc;
  _fSourcePosition_Y_Mpc = aSourcePosition.y()*invMpc;
  _fInitialPosition_Y_Mpc = aInitPosition.y()*invMpc;
  _fPosition_Y_Mpc = aPosition.y()*invMpc;
  _fSourcePosition_Z_Mpc = aSourcePosition.z()*invMpc;
  _fInitialPosition_Z_Mpc = aInitPosition.z()*invMpc;
  _fPosition_Z_Mpc = aPosition.z()*invMpc;
  _fEnergy_EeV = aMomentum.mag()*c_light*invEeV;
  _fMomentum_theta = aMomentum.theta();
  _fMomentum_phi = aMomentum.phi();
  _fSourceEnergy_EeV = aSourceEnergy*invEeV;
  for (int k=0;k<170;k++) (_fSpectrum)[k]=(float)(aSpectrum[k]);
  _fShowerTree->Fill();
}

void TRootOutput::Add1DNeutrino(int aType,
				double aSourcePosition,
				double aInitPosition,
				double aInitEnergy, 
				double aEnergy,
				double aSourceEnergy, int aInitType) {  
  _fNeutrinoNtuple->Fill(aType,  aInitType,
			 aSourcePosition * invMpc,
			 aInitPosition * invMpc,
			 aInitEnergy * invEeV,
			 aEnergy * invEeV,
			 aSourceEnergy * invEeV);
}

void TRootOutput::Add3DNeutrino(int aType,
				 TVector3D aSourcePosition,
				 TVector3D aInitPosition, 
				 TVector3D aPosition,
				 TVector3D aMomentum,
				 double aSourceEnergy, int aInitType) {
  _fNeutrinoNtuple->Fill(aType, aInitType,
			 aSourcePosition.x() * invMpc,
			 aSourcePosition.y() * invMpc,
			 aSourcePosition.z() * invMpc,
			 aInitPosition.x() * invMpc,
			 aInitPosition.y() * invMpc,
			 aInitPosition.z() * invMpc,
			 aPosition.x() * invMpc,
			 aPosition.y() * invMpc,
			 aPosition.z() * invMpc,
			 aMomentum.mag() * c_light * invEeV,
			 aMomentum.theta(),
			 aMomentum.phi(),
			 aSourceEnergy * invEeV);
}

TRootOutput::~TRootOutput() {
  _fRootFile->Write() ;
  _fRootFile->Close() ;
  cout << "ROOT file closed" << endl;
  //  delete _fEventNtuple ;
  //  delete _fRootFile ;
  // delete _fNeutrinoNtuple ;
}

//#endif

