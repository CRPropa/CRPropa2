
/**
   @file    lssfield.cc
   @author  Eric Armengaud, armengau@in2p3.fr
   @brief   Implementation of the TLSSField class. See the .h file
*/

#include "lssfield.h"

TLSSField::TLSSField(const char* aFileName) : TXmlParam(aFileName){
	SetType(MAGFIELD_LSS);
  // Grid parameters
  string lDum ;
  TiXmlElement* lpXmlField = XmlExtract().GetElement("MagneticField") ;
  TiXmlElement* lpXmlNx = lpXmlField->FirstChildElement("Nx") ;
  TiXmlElement* lpXmlNy = lpXmlField->FirstChildElement("Ny") ;
  TiXmlElement* lpXmlNz = lpXmlField->FirstChildElement("Nz") ;
  TiXmlElement* lpXmlStep = lpXmlField->FirstChildElement("Step_Mpc") ;
  if (!lpXmlNx || !lpXmlNy || !lpXmlNz || !lpXmlStep) {
    throw TCrpErr("Magnetic grid not properly specified in xml file"); 
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
      throw TCrpErr("Magnetic grid origin not properly specified");
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
  _fpB = new TVector3D[_fNx*_fNy*_fNz] ;
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
    lMagstream.getline(cha,100) ;
    lMagstream.getline(cha,100) ;
    _fBmax=0.;
    double lDumX, lDumY, lDumZ;
    for (int i=0; i<_fNx; i++) {
      for (int j=0; j<_fNy; j++) {
	for (int k=0; k<_fNz; k++) {
	  lMagstream >> lDumX >> lDumY >> lDumZ ;
	  (_fpB[i*_fNy*_fNz+j*_fNz+k]).set(lDumX,lDumY,lDumZ) ;
	  (_fpB[i*_fNy*_fNz+j*_fNz+k]) *= gauss ; // unit in the ascii file

	  long lInd = i*_fNy*_fNz+j*_fNz+k ;	  
	  //if( _fBmax < (_fpB[lInd]).mag() ){ std::cout<<"lInd"<<lInd<<"\t_fpB[lInd]).mag()="<<(_fpB[lInd]).mag()<<std::endl;}
	  _fBmax = max(_fBmax, (_fpB[lInd]).mag() ) ;
	}
      }
    }
    if (!lMagstream.good()) throw TCrpErr("Error reading the B field.") ;
    lMagstream.close() ;

  } else if (lFileType == "FITS") {

    TiXmlNode* lpXmlFileName = lpXmlFile->FirstChild() ;
    if ( !lpXmlFileName || !(lpXmlFileName->ToText()) ) {
      throw TXmlErr("Incorrect magnetic field file name");
    }
    string lFileName =lpXmlFileName->Value() ;
    int lStatus = 0 ;
    int lAnyNull = 0 ;
    int lNaxis = 0 ;
    long *lNaxes = new long[3] ; 

    double *lpBx = new double[_fNx*_fNy*_fNz] ; 
    fitsfile *lpBxFits ;
    //   const char *lExtBx = (lFileName+"[0]").c_str() ;
    // This does not seem to work anymore with cfitsio 3 :-( 
    const char *lExtBx = (lFileName).c_str() ;
    fits_open_file(&lpBxFits, lExtBx, READONLY, &lStatus) ;
    int lDummy = 0 ;
    fits_movabs_hdu(lpBxFits,1,&lDummy,&lStatus);
    fits_get_img_dim(lpBxFits, &lNaxis, &lStatus) ;
    if (lNaxis != 3) throw TFitsErr("NAXIS incorrect for Bx.") ;
    fits_get_img_size(lpBxFits, 3, lNaxes, &lStatus) ;
    if (lNaxes[0]!=_fNx || lNaxes[1]!=_fNy || lNaxes[2]!=_fNz) 
      throw TFitsErr("size of Bx array not ok with xml specifications.") ;
    fits_read_img(lpBxFits, TDOUBLE, 1, _fNx*_fNy*_fNz,
		  NULL, lpBx, &lAnyNull, &lStatus) ;
    fits_close_file(lpBxFits,&lStatus) ;
    if (lStatus) throw TFitsErr(lStatus) ;
    if (lAnyNull) throw TFitsErr("Undefined elements in Bx extension.") ;

    double *lpBy = new double[_fNx*_fNy*_fNz] ; 
    fitsfile *lpByFits ;
    //    const char *lExtBy = (lFileName+"[1]").c_str() ;
    const char *lExtBy = (lFileName).c_str() ;
    fits_open_file(&lpByFits, lExtBy, READONLY, &lStatus) ;
    fits_movabs_hdu(lpByFits,2,&lDummy,&lStatus);
    fits_get_img_dim(lpByFits, &lNaxis, &lStatus) ;
    if (lNaxis != 3) throw TFitsErr("NAXIS incorrect for By.") ;
    fits_get_img_size(lpByFits, 3, lNaxes, &lStatus) ;
    if (lNaxes[0]!=_fNx || lNaxes[1]!=_fNy || lNaxes[2]!=_fNz) 
      throw TFitsErr("size of By array not ok with xml specifications.") ;
    fits_read_img(lpByFits, TDOUBLE, 1, _fNx*_fNy*_fNz,
		  NULL, lpBy, &lAnyNull, &lStatus) ;
    fits_close_file(lpByFits,&lStatus) ;
    if (lStatus) throw TFitsErr(lStatus) ;
    if (lAnyNull) throw TFitsErr("Undefined elements in Bx extension.") ;

    double *lpBz = new double[_fNx*_fNy*_fNz] ; 
    fitsfile *lpBzFits ;
    const char *lExtBz = (lFileName).c_str() ;
    fits_open_file(&lpBzFits, lExtBz, READONLY, &lStatus) ;
    fits_movabs_hdu(lpBzFits,3,&lDummy,&lStatus);
    fits_get_img_dim(lpBzFits, &lNaxis, &lStatus) ;
    if (lNaxis != 3) throw TFitsErr("NAXIS incorrect for Bz.") ;
    fits_get_img_size(lpBzFits, 3, lNaxes, &lStatus) ;
    if (lNaxes[0]!=_fNx || lNaxes[1]!=_fNy || lNaxes[2]!=_fNz) 
      throw TFitsErr("size of Bx array not ok with xml specifications.") ;
    fits_read_img(lpBzFits, TDOUBLE, 1, _fNx*_fNy*_fNz,
		  NULL, lpBz, &lAnyNull, &lStatus) ;
    fits_close_file(lpBzFits,&lStatus) ;
    if (lStatus) throw TFitsErr(lStatus) ;
    if (lAnyNull) throw TFitsErr("Undefined elements in Bz extension.") ;

    _fBmax = 0. ;
    long lInd, lIndFits ;
    for (int i=0; i<_fNx; i++) {
      for (int j=0; j<_fNy; j++) {
	for (int k=0; k<_fNz; k++) {
	  lIndFits = i+j*_fNx+k*_fNx*_fNy ;
	  lInd = i*_fNy*_fNz+j*_fNz+k ;
	  (_fpB[lInd]).set(lpBx[lIndFits], lpBy[lIndFits], lpBz[lIndFits]);
	  (_fpB[lInd]) *= gauss ; // unit in the fits file
	  //if( _fBmax < (_fpB[lInd]).mag() ){ 
	  //std::cout<<"lInd"<<lInd<<"\t_fpB[lInd]).mag()="<<(_fpB[lInd]).mag()<<std::endl;
	  //}
	  _fBmax = max(_fBmax, (_fpB[lInd]).mag() ) ;
	}
      }
    }
    delete[] lpBx;
    delete[] lpBy;
    delete[] lpBz ;
    delete[] lNaxes ;

  } else throw TCrpErr("Error getting file type for magnetic field");
  
}

TLSSField::~TLSSField() {
  delete[] _fpB ;
}

