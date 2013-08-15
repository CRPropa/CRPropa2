

/**
   @file    sourcedensity.cc
   @author  Eric Armengaud, armengau@in2p3.fr
   @brief   Implementation of the TSourceDensity class. See the .h file
*/

#include "sourcedensity.h"

TSourceDensity::TSourceDensity(const char* aFileName) : TXmlParam(aFileName) {

  // Grid parameters
  string lDum ;

  TiXmlElement* lpXmlEnv = XmlExtract().GetElement("Environment") ;
  string lUnivType = lpXmlEnv->Attribute("type") ;
  TiXmlElement* lpXmlSources = XmlExtract().GetElement("Sources") ;
  TiXmlElement* lpXmlDensity = lpXmlSources->FirstChildElement("Density") ;
  string lDensityType = GetAttribute(lpXmlDensity,"type") ;

  if (lDensityType == "Uniform") { // Uniform density : no grid is needed
    _fIsUniform = 1 ;
    _fOrigin.set(0,0,0) ;
    _fNx = 0 ;
    _fNy = 0 ;
    _fNz = 0 ;
    _fStepsize = 0 ;
    _fpDensityArray = NULL ;
    _fpRandDistriX = NULL ;
    if (lUnivType == "One Dimension") { // get Xmax from the density
      TiXmlElement* lpXmlXmax = lpXmlDensity->FirstChildElement("Xmax_Mpc") ;
      if (!lpXmlXmax) throw TXmlErr("Xmax unspecified in uniform source density") ;
      lDum = lpXmlXmax->Attribute("value",&_fXmax) ;
      if(_fXmax > XmlExtract().GetDouble("MaxTime_Mpc") ) std::cout << "SourceDensity 1D larger than MaxTime_Mpc. Primaries may not reach the observer." << std::endl;
      _fXmax *= Mpc ;
      _fXmin = 0 ;
      TiXmlElement* lpXmlXmin = lpXmlDensity->FirstChildElement("Xmin_Mpc") ;
      // Xmin can be specified too
      if (lpXmlXmin) lDum = lpXmlXmin->Attribute("value",&_fXmin) ;
      _fXmin *= Mpc ;
      if (_fXmax <= _fXmin) throw TXmlErr("Source density: Xmax smaller than Xmin") ;
      _fYmin = 0 ;
      _fYmax = 0 ;
      _fZmin = 0 ;
      _fZmax = 0 ;
    } else { // get Xmax either from the density, or from the environment
      TiXmlElement* lpXmin = lpXmlEnv->FirstChildElement("Xmin_Mpc") ;
      TiXmlElement* lpBound = lpXmlEnv ;
      lpBound = lpXmlDensity ;
      lpXmin = lpBound->FirstChildElement("Xmin_Mpc") ;
      TiXmlElement* lpXmax = lpBound->FirstChildElement("Xmax_Mpc") ;
      TiXmlElement* lpYmin = lpBound->FirstChildElement("Ymin_Mpc") ;
      TiXmlElement* lpYmax = lpBound->FirstChildElement("Ymax_Mpc") ;
      TiXmlElement* lpZmin = lpBound->FirstChildElement("Zmin_Mpc") ;
      TiXmlElement* lpZmax = lpBound->FirstChildElement("Zmax_Mpc") ;
      if (!lpXmin || !lpXmax || !lpYmin || !lpYmax || !lpZmin || !lpZmax)
	throw TXmlErr("Density boundary missing.") ;
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

  } else if (lDensityType == "Grid") {
    _fIsUniform = 0 ;
    TiXmlElement* lpXmlNx = lpXmlDensity->FirstChildElement("Nx") ;
    TiXmlElement* lpXmlNy = lpXmlDensity->FirstChildElement("Ny") ;
    TiXmlElement* lpXmlNz = lpXmlDensity->FirstChildElement("Nz") ;
    TiXmlElement* lpXmlStep = lpXmlDensity->FirstChildElement("Step_Mpc") ;
    
    if (!lpXmlNx || ( (!lpXmlNy || !lpXmlNz) && lUnivType != "One Dimension" ) 
	|| !lpXmlStep) throw TCrpErr("Source density not properly specified in xml file"); 
    lDum = lpXmlNx->Attribute("value",&_fNx) ;
    
    if (lUnivType != "One Dimension") {
      lDum = lpXmlNy->Attribute("value",&_fNy) ;
      lDum = lpXmlNz->Attribute("value",&_fNz) ;
    } else {
      _fNy = 1 ; // 1D array in this case
      _fNz = 1 ;
    }
    lDum = lpXmlStep->Attribute("value",&_fStepsize) ;
    _fStepsize *= Mpc ;
    
    TiXmlElement* lpXmlOrigin = lpXmlDensity->FirstChildElement("Origin") ;
    if (lpXmlOrigin) {
      TiXmlElement* lpXmlOx = lpXmlOrigin->FirstChildElement("X_Mpc") ;
      TiXmlElement* lpXmlOy = lpXmlOrigin->FirstChildElement("Y_Mpc") ;
      TiXmlElement* lpXmlOz = lpXmlOrigin->FirstChildElement("Z_Mpc") ;
      if (!lpXmlOx || !lpXmlOy || !lpXmlOz) {
	throw TCrpErr("Density grid origin not properly specified");
      }
      double lOx,lOy,lOz;
      lDum = lpXmlOx->Attribute("value",&lOx) ;
      lDum = lpXmlOy->Attribute("value",&lOy) ;
      lDum = lpXmlOz->Attribute("value",&lOz) ;
      _fOrigin.set(lOx*Mpc,lOy*Mpc,lOz*Mpc) ;
    } else {
      _fOrigin.set(0,0,0) ;
    }

    _fXmin = _fOrigin.x() ;
    _fXmax = _fXmin + _fNx * _fStepsize ;
    _fYmin = _fOrigin.y() ;
    _fYmax = _fYmin + _fNy * _fStepsize ;
    _fZmin = _fOrigin.z() ;
    _fZmax = _fZmin + _fNz * _fStepsize ;

    // Filling the array
    _fpDensityArray = new double[_fNx*_fNy*_fNz] ;
    TiXmlElement* lpXmlFile = lpXmlDensity->FirstChildElement("File") ;
    string lFileType = lpXmlFile->Attribute("type") ;
    TiXmlNode* lpXmlFileName = lpXmlFile->FirstChild() ;
    if ( !lpXmlFileName || !(lpXmlFileName->ToText()) )
      throw TCrpErr("Error parsing xml : source density grid file name");
    string lFileName =lpXmlFileName->Value() ;
    
    if (lFileType == "ASCII") {
      fstream lDensitystream(lFileName.c_str(), ios::in) ;
      if (! lDensitystream.is_open())
	throw TCrpErr("Error opening density grid file : " + lFileName );
      char cha[256] ;
      lDensitystream.getline(cha,100) ;
      lDensitystream.getline(cha,100) ;
      double lDumRho;
      for (int i=0; i<_fNx; i++) {
	for (int j=0; j<_fNy; j++) {
	  for (int k=0; k<_fNz; k++) {
	    lDensitystream >> lDumRho;
	    (_fpDensityArray[i*_fNy*_fNz+j*_fNz+k]) = lDumRho ; // density is not normalized
	  }
	}
      }
      if (!lDensitystream.good()) throw TCrpErr("Error reading the density grid file.") ;
      lDensitystream.close() ;
      
      
    } else if (lFileType == "FITS") {
      int lStatus = 0 ;
      int lAnyNull = 0 ;
      int lNaxis = 0 ;
      double *lpDensity = new double[_fNx*_fNy*_fNz] ; 
      fitsfile *lpDensityFits ;
      string lStrTmp = lFileName ;
      const char *lExtDensity = lStrTmp.c_str() ;
      fits_open_file(&lpDensityFits, lExtDensity, READONLY, &lStatus) ;
      fits_get_img_dim(lpDensityFits, &lNaxis, &lStatus) ;
      if (lUnivType != "One Dimension" && lNaxis != 3)
	throw TFitsErr("NAXIS incorrect for 3D density.") ;
      if (lUnivType == "One Dimension" && lNaxis != 1) 
	throw TFitsErr("NAXIS incorrect for 1D density.") ;
      if (lUnivType != "One Dimension") {
	long *lNaxes = new long[3] ; 
	fits_get_img_size(lpDensityFits, 3, lNaxes, &lStatus) ;
	if (lNaxes[0]!=_fNx || lNaxes[1]!=_fNy || lNaxes[2]!=_fNz) 
	  throw TFitsErr("size of density array not ok with xml specifications.") ;
	delete[] lNaxes ;
      } else {
	long *lNaxes = new long[1] ; 
	fits_get_img_size(lpDensityFits, 1, lNaxes, &lStatus) ;
	if (lNaxes[0]!=_fNx) 
	  throw TFitsErr("size of density array not ok with xml specifications.") ;
	delete[] lNaxes ;
      }
      fits_read_img(lpDensityFits, TDOUBLE, 1, _fNx*_fNy*_fNz,
		    NULL, lpDensity, &lAnyNull, &lStatus) ;
      fits_close_file(lpDensityFits,&lStatus) ;
      if (lStatus) throw TFitsErr(lStatus) ;
      if (lAnyNull) throw TFitsErr("Undefined elements in density extension.") ;
      
      long lInd, lIndFits ;
      for (int i=0; i<_fNx; i++) {
	for (int j=0; j<_fNy; j++) {
	  for (int k=0; k<_fNz; k++) {
	    lIndFits = i+j*_fNx+k*_fNx*_fNy ;
	    lInd = i*_fNy*_fNz+j*_fNz+k ; // long() pour etre sur que c'est bon!
	    _fpDensityArray[lInd] = lpDensity[lIndFits] ; // no normalization
	  }
	}
      }
      delete[] lpDensity;
      
    } else throw TCrpErr("Error getting file type for source density grid" );
    
    // RandGeneral object associated with the density
    // 1) In the case of 1D, the density array has to be binned enough, 
    // so that sources will be drawn uniformly in each "bin" of the array
    // 2) In 3D, we iteratively select each of the 3 indices, and then build a subarray...
    int lIntType = 1 ;
    if (lUnivType == "One Dimension") {
      lIntType = 0 ;
      _fpRandDistriX = new RandGeneral(_fpDensityArray,_fNx,lIntType) ;
      // The seed for this object is the same as for the RandFlat arrays (ie. can be specified)
    } else {
      double lDdum ;
      double *lpDensityX = new double[_fNx] ;
      double *lpDensityXY = new double[_fNy] ;
      for (int i=0; i<_fNx; i++) {
	lpDensityX[i] = 0 ;
	for (int j=0; j<_fNy; j++) {
	  lpDensityXY[j] = 0 ;
	  for (int k=0; k<_fNz; k++) {
	    lDdum = _fpDensityArray[i*_fNy*_fNz+j*_fNz+k] ;
	    lpDensityX[i] += lDdum ;
	    lpDensityXY[j] += lDdum ;
	  }
	}
	_fpRandDistriXY.push_back(new RandGeneral(lpDensityXY,_fNy,lIntType)) ;
      }
      _fpRandDistriX = new RandGeneral(lpDensityX,_fNx,lIntType) ;
      delete[] lpDensityX ;
      delete[] lpDensityXY ;
    }

  } else throw TXmlErr("Unknown source density type") ;

}


TSourceDensity::~TSourceDensity() {

  delete[] _fpDensityArray ;
  delete _fpRandDistriX ;
  while (_fpRandDistriXY.size()) {
    delete _fpRandDistriXY.back() ;
    _fpRandDistriXY.pop_back() ;
  }

}

double TSourceDensity::getDensity(TVector3D aPosition) const {
  // Subroutine to interpolate the density at any position
  // aPosition must be inside a shell around the grid

  if (_fIsUniform) throw TCrpErr("Attempting to evaluate the source density without density grid.") ;

  double lLocDensity ;
  int lQx,lQy,lQz ;
  int lQx1,lQy1,lQz1 ;
  double a,b,c,d,e,f ;
  TVector3D lPosInStep = (aPosition-_fOrigin)/_fStepsize ;
  lQx = (int)floor(lPosInStep.x()) ;
  a=lPosInStep.x()-lQx ;
  d=1-a ;
  lQx1 = lQx+1 ;
  lQx = lQx % _fNx ;
  if (lQx<0) lQx += _fNx ;
  lQx1 = lQx1 % _fNx ;
  if (lQx1<0) lQx1 += _fNx ;
  lQy = (int)floor(lPosInStep.y()) ;
  b=lPosInStep.y()-lQy ;
  e=1-b ;
  lQy1 = lQy+1 ;
  lQy = lQy % _fNy ;
  if (lQy<0) lQy += _fNy ;
  lQy1 = lQy1 % _fNy ;
  if (lQy1<0) lQy1 += _fNy ;
  lQz = (int)floor(lPosInStep.z()) ;
  c=lPosInStep.z()-lQz ;
  f=1-c ;
  lQz1 = lQz+1 ;
  lQz = lQz % _fNz ;
  if (lQz<0) lQz += _fNz ;
  lQz1 = lQz1 % _fNz ;
  if (lQz1<0) lQz1 += _fNz ;

  int NyNz = _fNy*_fNz; 
  lLocDensity = 
    d*(e*(f*_fpDensityArray[lQx*NyNz+lQy*_fNz+lQz]+c*_fpDensityArray[lQx*NyNz+lQy*_fNz+lQz1])
       +b*(f*_fpDensityArray[lQx*NyNz+lQy1*_fNz+lQz]+c*_fpDensityArray[lQx*NyNz+lQy1*_fNz+lQz1]))+
    a*(e*(f*_fpDensityArray[lQx1*NyNz+lQy*_fNz+lQz]+c*_fpDensityArray[lQx1*NyNz+lQy*_fNz+lQz1])+
       b*(f*_fpDensityArray[lQx1*NyNz+lQy1*_fNz+lQz]+c*_fpDensityArray[lQx1*NyNz+lQy1*_fNz+lQz1])) ;

  return lLocDensity;
}

TVector3D TSourceDensity::getSourcePosition() const {

  double lPosX ;
  TVector3D lPosition ;
  string lUnivType = (XmlExtract().GetElement("Environment"))->Attribute("type") ;
  if (lUnivType == "One Dimension") {

    if (this->IsUniform()) {
      lPosX = _fXmin + RandFlat::shoot()*(_fXmax-_fXmin) ;
    } else {
      lPosX = _fpRandDistriX->shoot() ;
      lPosX *= (_fNx-1)*_fStepsize ;
      lPosX += _fOrigin.x() ;
    }
    lPosition.set(lPosX,0,0) ;

  } else {

    if (this->IsUniform()) {
      lPosition.set(_fXmin+(_fXmax-_fXmin)*RandFlat::shoot(),
		    _fYmin+(_fYmax-_fYmin)*RandFlat::shoot(),
		    _fZmin+(_fZmax-_fZmin)*RandFlat::shoot()) ;
    } else {
      // 1) Get indices (i,j,k) of the cell
      int lI,lJ,lK ;
      TVector3D lPos0, lDum ;
      double *lpDensityZ = new double[_fNz] ;
      lI = int(_fpRandDistriX->shoot()*_fNx) ;
      lJ = int((_fpRandDistriXY[lI])->shoot()*_fNy) ;
      for (int k=0; k<_fNz; k++) lpDensityZ[k] = _fpDensityArray[lI*_fNy*_fNz+lJ*_fNz+k] ;
      RandGeneral lRandDistriZ(lpDensityZ,_fNz,1) ;
      delete[] lpDensityZ ;
      lK = int(lRandDistriZ.shoot()*_fNz) ;
      lDum.set(lI,lJ,lK) ;
      lPos0 = _fOrigin + _fStepsize*lDum ;
      
      // 2) Choose a position inside the cell using a flat distribution
      lDum.set(RandFlat::shoot()-0.5,RandFlat::shoot()-0.5,RandFlat::shoot()-0.5) ;
      lPosition = lPos0 + _fStepsize * lDum ;
      
      // Check that position is well inside the box
      if (lPosition.x() <= _fXmin) lPosition.setX(lPosition.x()+_fStepsize*_fNx) ;
      if (lPosition.x() >= _fXmax) lPosition.setX(lPosition.x()-_fStepsize*_fNx) ;
      if (lPosition.y() <= _fYmin) lPosition.setY(lPosition.y()+_fStepsize*_fNy) ;
      if (lPosition.y() >= _fYmax) lPosition.setY(lPosition.y()-_fStepsize*_fNy) ;
      if (lPosition.z() <= _fZmin) lPosition.setZ(lPosition.z()+_fStepsize*_fNz) ;
      if (lPosition.z() >= _fZmax) lPosition.setZ(lPosition.z()-_fStepsize*_fNz) ;
    }

  }

  return lPosition ;
}

