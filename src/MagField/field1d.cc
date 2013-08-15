
/**
   @file    field1d.cc
   @author  Eric Armengaud, armengau@in2p3.fr
   @brief   Implementation of the TField1D class. See the .h file
*/

#include "field1d.h"

TField1D::TField1D(const char* aFileName) : TXmlParam(aFileName) {
  SetType(MAGFIELD_1D);
  TiXmlElement* lpXmlField = XmlExtract().GetElement("MagneticField") ;
  TiXmlElement* lpXmlFile = lpXmlField->FirstChildElement("File") ;
  string lFileType = lpXmlFile->Attribute("type") ;
  if (lFileType == "ASCII") {

    TiXmlNode* lpXmlFileName = lpXmlFile->FirstChild() ;
    if ( !lpXmlFileName || !(lpXmlFileName->ToText()) ) {
      throw TXmlErr("Incorrect magnetic field file name");
    }
    string lFileName =lpXmlFileName->Value() ;
    fstream lMagstream(lFileName.c_str(), ios::in) ;
    if (! lMagstream.is_open()) {
      throw TCrpErr("Error opening magnetic field file : " + lFileName );
    }
    char cha[256] ;
    lMagstream.getline(cha,100) ; // 2 lines of comment
    lMagstream.getline(cha,100) ;
    double lDumX, lDumB; // convention: X in Mpc, B in muG ! 
    while (lMagstream >> lDumX >> lDumB) {
      _fPositions.push_back(lDumX * Mpc) ;
      _fFieldValues.push_back(lDumB * muG) ;
      if (!lMagstream.good()) throw TCrpErr("Error reading magnetic file.") ;
    }
    lMagstream.close() ;

  } else throw TXmlErr("Error getting file type for magnetic field" );

  _fXmax = _fPositions[_fPositions.size()-1] ;
  _fStepsize = _fPositions[1] - _fPositions[0] ;
  if (_fStepsize <= 0) 
    throw TCrpErr("Stepsize in B grid is null or negative.") ;
  for (unsigned int i=0; i<_fPositions.size()-1; i++) {
    if ( abs(_fPositions[i+1]-_fPositions[i]-_fStepsize)/_fStepsize >= CHECK_STEP ) throw TCrpErr("Stepsize not constant in B grid.") ;
  }

}

TVector3D TField1D::getField(TVector3D aPosition) const {
  // as the accuracy of magnetic field topology is of no importance,
  // we just make a linear interpolation of the arrays.

  if (aPosition.x() < _fPositions[0] || aPosition.x() >= _fXmax) 
    throw TCrpErr("Trying to interpolate B field out of its bounds in TField1D") ;

  unsigned int i0 = 0 ;
  while (_fPositions[i0+1] <= aPosition.x() ) i0 += 1 ;
  double lFieldValue = _fFieldValues[i0] + (_fFieldValues[i0+1]-_fFieldValues[i0])
    *(aPosition.x()-_fPositions[i0])/_fStepsize ;

  TVector3D lField(0,lFieldValue,0) ; // convention : champ selon (Oy)
  return lField ;
}
