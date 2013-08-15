/**
   @file    lssgas.cc
   @author  Luca Maccione, luca.maccione@desy.de
   @brief   Implementation of the TLSSGas class. See the .h file
*/

#include "lssgas.h"

TLSSGas::TLSSGas(const char* aFileName) : TGridGas(), TXmlParam(aFileName){
  SetType(LSSGAS);
  // Grid parameters
  string lDum ;
  TiXmlElement* lpXmlField = XmlExtract().GetElement("Gas") ;
  TiXmlElement* lpXmlNx = lpXmlField->FirstChildElement("Nx") ;
  TiXmlElement* lpXmlNy = lpXmlField->FirstChildElement("Ny") ;
  TiXmlElement* lpXmlNz = lpXmlField->FirstChildElement("Nz") ;
  TiXmlElement* lpXmlStep = lpXmlField->FirstChildElement("Step_Mpc") ;
  if (!lpXmlNx || !lpXmlNy || !lpXmlNz || !lpXmlStep) {
    throw TCrpErr("Gas grid not properly specified in xml file"); 
  }
  lDum = lpXmlNx->Attribute("value",&_fNx) ;
  lDum = lpXmlNy->Attribute("value",&_fNy) ;
  lDum = lpXmlNz->Attribute("value",&_fNz) ;
  lDum = lpXmlStep->Attribute("value",&_fStepsize) ;
  _fStepsize *= Mpc ;
  TiXmlElement* lpXmlOrigin = lpXmlField->FirstChildElement("Origin") ;
  if (lpXmlOrigin) {
    TiXmlElement* lpXmlOx = lpXmlOrigin->FirstChildElement("X_Mpc") ;
    TiXmlElement* lpXmlOy = lpXmlOrigin->FirstChildElement("Y_Mpc") ;
    TiXmlElement* lpXmlOz = lpXmlOrigin->FirstChildElement("Z_Mpc") ;
    if (!lpXmlOx || !lpXmlOy || !lpXmlOz) {
      throw TCrpErr("Gas grid origin not properly specified");
    }
    double lOx,lOy,lOz;
    lDum = lpXmlOx->Attribute("value",&lOx) ;
    lDum = lpXmlOy->Attribute("value",&lOy) ;
    lDum = lpXmlOz->Attribute("value",&lOz) ;
    _fOrigin.set(lOx*Mpc,lOy*Mpc,lOz*Mpc) ;
  } else {
    _fOrigin = TVector3D() ;
  }

  // Filling the array
  _fpGas = vector<double>(_fNx*_fNy*_fNz, 0.0);
  TiXmlElement* lpXmlFile = lpXmlField->FirstChildElement("File") ;
  string lFileType = lpXmlFile->Attribute("type") ;
  if (lFileType == "ASCII") {

    TiXmlNode* lpXmlFileName = lpXmlFile->FirstChild() ;
    if ( !lpXmlFileName || !(lpXmlFileName->ToText()) ) {
      throw TXmlErr("Incorrect gas field file name");
    }
    string lFileName =lpXmlFileName->Value() ;
    fstream lGasstream(lFileName.c_str(), ios::in) ;
    if (! lGasstream.is_open()) {
      throw TCrpErr("Error opening gas field file : " + lFileName );
    }
    char cha[256] ;
    lGasstream.getline(cha,100) ;
    lGasstream.getline(cha,100) ;
    double lDumX;
    for (int i=0; i<_fNx; i++) {
      for (int j=0; j<_fNy; j++) {
	for (int k=0; k<_fNz; k++) {
	  lGasstream >> lDumX;
	  _fpGas[i*_fNy*_fNz+j*_fNz+k] = lDumX ;  // cm^-3
	}
      }
    }
    if (!lGasstream.good()) throw TCrpErr("Error reading the Gas field.") ;
    lGasstream.close() ;

  } else if (lFileType == "FITS") {

    TiXmlNode* lpXmlFileName = lpXmlFile->FirstChild() ;
    if ( !lpXmlFileName || !(lpXmlFileName->ToText()) ) {
      throw TXmlErr("Incorrect Gas field file name");
    }
    string lFileName =lpXmlFileName->Value() ;
    int lStatus = 0 ;
    int lAnyNull = 0 ;
    int lNaxis = 0 ;
    long *lNaxes = new long[3] ; 

    double *lpGas = new double[_fNx*_fNy*_fNz] ; 
    fitsfile *lpGasFits ;
    //   const char *lExtBx = (lFileName+"[0]").c_str() ;
    // This does not seem to work anymore with cfitsio 3 :-( 
    const char *lExtGas = (lFileName).c_str() ;
    fits_open_file(&lpGasFits, lExtGas, READONLY, &lStatus) ;
    int lDummy = 0 ;
    fits_movabs_hdu(lpGasFits,1,&lDummy,&lStatus);
    fits_get_img_dim(lpGasFits, &lNaxis, &lStatus) ;
    if (lNaxis != 3) throw TFitsErr("NAXIS incorrect for Gas.") ;
    fits_get_img_size(lpGasFits, 3, lNaxes, &lStatus) ;
    if (lNaxes[0]!=_fNx || lNaxes[1]!=_fNy || lNaxes[2]!=_fNz) 
      throw TFitsErr("size of Gas array not ok with xml specifications.") ;
    fits_read_img(lpGasFits, TDOUBLE, 1, _fNx*_fNy*_fNz,
		  NULL, lpGas, &lAnyNull, &lStatus) ;
    fits_close_file(lpGasFits,&lStatus) ;
    if (lStatus) throw TFitsErr(lStatus) ;
    if (lAnyNull) throw TFitsErr("Undefined elements in Gas extension.") ;

    _fGasmax = 0 ;
    long lInd, lIndFits ;
    for (int i=0; i<_fNx; i++) {
      for (int j=0; j<_fNy; j++) {
	for (int k=0; k<_fNz; k++) {
	  lIndFits = i+j*_fNx+k*_fNx*_fNy ;
	  lInd = i*_fNy*_fNz+j*_fNz+k ;
	  _fpGas[lInd] = lpGas[lIndFits];
	  _fGasmax = max(_fGasmax, _fpGas[lInd]) ;
	}
      }
    }
    delete[] lpGas;
    delete[] lNaxes ;

  } else throw TCrpErr("Error getting file type for Gas field");

}
