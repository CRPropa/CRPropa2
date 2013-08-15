/**
   @file    variableinfrared.cc
   @author  Luca Maccione, luca.maccione@desy.de
   @brief   Implementation of the TVariableInfrared class. See the .h file
*/

#include "variableinfrared.h"

using namespace std;

TVariableInfrared::TVariableInfrared(const char* aFileName, bool oned_) : IR(), TXmlParam(aFileName) {
  oned=oned_;
  string lDum ;
  string lSpectralShape = "Shells" ;
  bool lCheckIR = XmlExtract().CheckElement("InfraredBackground") ;
  TiXmlElement* lpXmlIR = NULL;
  if (lCheckIR) {
    lpXmlIR = XmlExtract().GetElement("InfraredBackground") ;
    if (lpXmlIR) {
      TiXmlElement* lpXmlSpectrum = lpXmlIR->FirstChildElement("SpectralShape") ;
      if (lpXmlSpectrum) lSpectralShape = lpXmlSpectrum->Attribute("type") ;
      TiXmlElement* lpXmlZmax = lpXmlIR->FirstChildElement("MaxRedshift") ;
      if (lpXmlZmax) lDum = lpXmlZmax->Attribute("value",&_fZmax) ; 
      
      
      if (lSpectralShape == "Shells") {
	SetType(SHELL);
      } else if (lSpectralShape == "3D") SetType(THREED);
      else throw TXmlErr("Unknown IR background type") ;
      
      TiXmlElement* lpXmlFile = lpXmlIR->FirstChildElement("File");
      string lFileType = lpXmlFile->Attribute("type");
      TiXmlNode* lpXmlFileName = lpXmlFile->FirstChild();
      if ( !lpXmlFileName || !(lpXmlFileName->ToText()) )
	throw TCrpErr("Error parsing xml : IR density grid file name");
      string lFileName = lpXmlFileName->Value();
      //      cout << lFileName << endl;
      if (Type() == SHELL) {
	if (lFileType == "ASCII") {
	  fstream lDensitystream(lFileName.c_str(), ios::in);
	  if (!lDensitystream.is_open())
	    throw TCrpErr("Error opening IR grid file : " + lFileName);
	  
	  string str;
	  getline(lDensitystream, str);
	  
	  readstring(str, _fPositionX);
	  _fEnergyDensity = vector< vector<double> >(_fPositionX.size(), vector<double>());
	  _fNx = _fPositionX.size();
	  while (lDensitystream.good()) {
	    str.clear();  // Just to be sure
	    getline(lDensitystream,str);
	    readstring(str, _fEnergy, _fEnergyDensity);
	  }
	  _fEnergy.pop_back();
	  _fNE = _fEnergy.size();
	  
	  lDensitystream.close();
	  
	  for (int i = 0; i < _fNx-1; i++) {
	    for (int j = 0; j < _fNE; j++) {
	      _fEnergyDensity[i][j] += _fEnergyDensity[_fNx-1][j];
// 	      static int testflag=0;
// 	      if(!testflag) {std::cout<< "variable Infrared modified!!!"<< std::endl; testflag=1;}
// 	      _fEnergyDensity[i][j] = _fEnergyDensity[_fNx-1][j];
// 	      cout << "EBL = " << _fEnergyDensity[_fNx-1][j] << endl; // For debugging
	    }
	  }
	  //	cout << _fEnergy.front() << " " << _fEnergy.back() << " " << _fNE << endl;
	  if (_fEnergy.front() > _fEnergy.back()) {
	    
	    reverse(_fEnergy.begin(), _fEnergy.end());
	    for (int i = 0; i < _fNx; i++) reverse(_fEnergyDensity[i].begin(), _fEnergyDensity[i].end());
	  }
	  
	}
	else if (lFileType == "FITS") {
	  throw TCrpErr("FITS files in SHELL IR spectral type not implemented yet");
	}
      }
      else if (Type() == THREED) {
	if (lFileType == "ASCII") {
	  fstream lDensitystream(lFileName.c_str(), ios::in);
	  if (!lDensitystream.is_open())
	    throw TCrpErr("Error opening IR grid file : " + lFileName);

	  string str;
	  getline(lDensitystream, str);
	  readstring(str, _fPositionX);
	  _fNx = _fPositionX.size();
	  str.clear();

	  getline(lDensitystream, str);
	  readstring(str, _fPositionY);
	  _fNy = _fPositionY.size();
	  str.clear();

	  getline(lDensitystream, str);
	  readstring(str, _fPositionZ);
	  _fNz = _fPositionZ.size();
	  str.clear();

	  getline(lDensitystream, str);
	  readstring(str, _fEnergy);
	  _fNE = _fEnergy.size();
	  str.clear();
	  
	  _fEnergyDensity = vector< vector<double> >(_fNx*_fNy*_fNz, vector<double>(_fNE,0.0));
	  int counterPos = 0;
	  int counterE = 0;
	  while (lDensitystream.good() && counterE*counterPos < _fNx*_fNy*_fNz*_fNE) {
	    lDensitystream >> _fEnergyDensity[counterPos][counterE];
	    counterE++;
	    if (counterE%_fNE == 0) {
	      counterE = 0;
	      counterPos++;
	    }
	  }
	  
	  lDensitystream.close();
	  
	}
	else if (lFileType == "FITS") {
	  throw TCrpErr("FITS files in THREED IR spectral type not implemented yet");
	}
      }
      else throw TCrpErr("Unknown variable photon field spectral type.");
    }
  }
  else throw TCrpErr("Error getting file type for IR background");
}

TPhotonSpectrum TVariableInfrared::Spectrum(double x, double y, double z) {

  if (Type() == SHELL) {
    
    if (x >= _fPositionX[_fNx-2]) return TPhotonSpectrum(_fEnergy, vector<double>(_fEnergyDensity[_fNx-1]));
    else if (x <= _fPositionX.front()) return TPhotonSpectrum(_fEnergy, vector<double>(_fEnergyDensity[0]));
    
    
    // Assume linear spacing in r
    int j = int( floor( (x-_fPositionX.front())/(_fPositionX[1]-_fPositionX.front()) ) );
    if (j == _fNx-1) j--;

    int i = j+1;
    double r1 = (x-_fPositionX[j])/(_fPositionX[i]-_fPositionX[j]);

    vector<double> result(_fNE, 0.0);
    for (int k = 0; k < _fNE; k++) result[k] = (1.0-r1)*_fEnergyDensity[j][k] + r1 * _fEnergyDensity[i][k];
    
    return TPhotonSpectrum(_fEnergy, result);
  }
  else {
    int i = int( floor( (x-_fPositionX.front())/(_fPositionX[1]-_fPositionX.front()) ) );
    int j = int( floor( (y-_fPositionY.front())/(_fPositionY[1]-_fPositionY.front()) ) );
    int k = int( floor( (z-_fPositionZ.front())/(_fPositionZ[1]-_fPositionZ.front()) ) );

    double xd = (x-_fPositionX[i])/(_fPositionX[i+1]-_fPositionX[i]);
    double yd = (y-_fPositionY[j])/(_fPositionY[j+1]-_fPositionY[j]);
    double zd = (z-_fPositionZ[k])/(_fPositionZ[k+1]-_fPositionZ[k]);
    
    vector<double> result(_fNE, 0.0);
    for (int ee = 0; ee < _fNE; ee++) {
      double i1 = _fEnergyDensity[(_fNy*i+j)*_fNz+k][ee]*(1.0-zd) + zd*_fEnergyDensity[(_fNy*i+j)*_fNz+k+1][ee];
      double i2 = _fEnergyDensity[(_fNy*i+j+1)*_fNz+k][ee]*(1.0-zd) + zd*_fEnergyDensity[(_fNy*i+j+1)*_fNz+k+1][ee];
      double j1 = _fEnergyDensity[(_fNy*(i+1)+j)*_fNz+k][ee]*(1.0-zd) + zd*_fEnergyDensity[(_fNy*(i+1)+j)*_fNz+k+1][ee];
      double j2 = _fEnergyDensity[(_fNy*(i+1)+j+1)*_fNz+k][ee]*(1.0-zd) + zd*_fEnergyDensity[(_fNy*(i+1)+j+1)*_fNz+k+1][ee];
      
      double w1 = i1*(1.0-yd)+i2*yd;
      double w2 = j1*(1.0-yd)+j2*yd;
      
      result[ee] = w1*(1.0-xd)+xd*w2;
    }
    return TPhotonSpectrum(_fEnergy, result);
  }

}

double TVariableInfrared::GetPhotonDensity(double x, 
					   double y, 
					   double z, 
					   double redshift,
					   double PhotonEnergy /* GeV */) {

  TPhotonSpectrum sp = Spectrum(x,y,z);
  PhotonEnergy=PhotonEnergy*1.e9;
  /*
    The units of TPhotonSpectrum are eV/cm^3 = 1.e-9*1.e6 GeV/m^3. The is because the density is given as E^2 * dN/dE/dV. Needed are the same units as in the CMB case: #/Gev/m^3 - because then integration over the  energy gives number per m^3.  Thus, division by (E/GeV)^2 is needed. After that one has to change the unit system to h_bar=c=1. In this units one should end up with GeV^2. Thus the result has to be multiplied with a factor of (h_bar^3*c)^3 which has units of (GeV*m)^3. */
  
  // return value in GeV^2:                           ( h_bar     * c    )^3      *  1.e-9 * 1.e6 / pow(1e-9*PhotonEnergy,2)
  return (
	  sp.Spectrum(PhotonEnergy) /*eV/cm^3 */ * pow( 6.528e-25 * 3.e8 , 3 )    *  1.e-9 * 1.e6 / pow(1e-9*PhotonEnergy,2)
	  ); 
}

void readstring(const string& str, vector<double>& en, vector< vector<double> >& spectrum) {

  size_t pos = 0;
  pos = str.find(" ", pos);
  size_t pos1 = str.find(" ", pos+1);
  while (pos1 == pos+1) {
    pos = pos1;
    pos1 = str.find(" ", pos+1);
  }
  string value(str.substr(pos+1,pos1-1-pos));
  en.push_back(atof(value.c_str()));
  //  cout << "Energy = " << en.back() << endl;
  pos = pos1-1;
 
  value.clear();
  int counter = 0;

  while (pos < str.size()) {
    pos = str.find(" ", pos);
    pos1 = str.find(" ", pos+1);
    while (pos1 == pos+1) {
      pos = pos1;
      pos1 = str.find(" ", pos+1);
    }
    value = string(str.substr(pos+1,pos1-1-pos));
    spectrum[counter].push_back(atof(value.c_str()));
    //   cout << "Spectrum = " << spectrum[counter].back() << endl;
    pos = pos1-1;
    counter++;
  }

  return ;
}

void readstring(const string& str, vector<double>& r) {

  size_t pos = 0;
  pos = str.find(" ", pos);
  size_t pos1 = str.find(" ", pos+1);
  while (pos1 == pos+1) {
    pos = pos1;
    pos1 = str.find(" ", pos+1);
  }
  string value(str.substr(pos+1,pos1-1-pos));
  pos = pos1-1;

  while(pos < str.size()) {
    pos = str.find(" ", pos);
    pos1 = str.find(" ", pos+1);
    while (pos1 == pos+1) {
      pos = pos1;
      pos1 = str.find(" ", pos+1);
    }
    string value(str.substr(pos+1,pos1-1-pos));
    r.push_back(atof(value.c_str()));
    //    cout << r.back() << endl;
    pos = pos1-1;
  }

  return ;
}
