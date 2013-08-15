/**
   @file    clustergas.cc
   @author  Luca Maccione, luca.maccione@desy.de
   @brief   Implementation of the TCLusterGas class. See the .h file
*/

#include "clustergas.h"

using namespace std;

TClusterGas::TClusterGas(const char* aFileName) : TXmlParam(aFileName) {
  SetType(CLUSTER);

  bool lCheckGas = XmlExtract().CheckElement("Gas") ;
  TiXmlElement* lpXmlGas = NULL;
  if (lCheckGas) {
    lpXmlGas = XmlExtract().GetElement("Gas") ;
    TiXmlElement* lpXmlFile = lpXmlGas->FirstChildElement("File");
    string lFileType = lpXmlFile->Attribute("type");
    TiXmlNode* lpXmlFileName = lpXmlFile->FirstChild();
    if ( !lpXmlFileName || !(lpXmlFileName->ToText()) )
      throw TCrpErr("Error parsing xml : Gas density grid file name");
    string lFileName = lpXmlFileName->Value();
 
    if (lFileType == "ASCII") {
      fstream lDensitystream(lFileName.c_str(), ios::in);
      if (!lDensitystream.is_open())
	throw TCrpErr("Error opening Gas grid file : " + lFileName);
      
      double r = 0.0;
      double dens = 0.0;
      while (lDensitystream.good()) {
	lDensitystream >> r >> dens;
	_fRarray.push_back(r);
	_fDensity.push_back(dens);
      }
      _fRarray.pop_back();
      _fDensity.pop_back();
      lDensitystream.close();
    }
  } else throw TCrpErr("Error getting file type for Gas background"); 
}

double TClusterGas::Density(double r) {
  if (r > _fRarray.back()) return _fDensity.back();

  int i = 0;
  while (_fRarray[i] < r) i++;
  int j = i-1;
  double r1 = (r-_fRarray[j])/(_fRarray[i]-_fRarray[j]);

  double result = (1.0-r1)*_fDensity[j] + r1 * _fDensity[i];

  return result;
}
