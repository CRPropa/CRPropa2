
/**
   @file    basicpinteractions.cc
   @author  Eric Armengaud, armengau@in2p3.fr
   @brief   Implementation of the TBasicPInteractions class. See the .h file
*/

#include "basicpinteractions.h"

TBasicPInteractions::TBasicPInteractions(const char* aFileName) : TXmlParam(aFileName) {
SetType(INTERACTION_BASIC);
  string lIntDir;

  TiXmlElement* lpXmlInter = XmlExtract().GetElement("Interactions") ;

  // Interaction directory
  TiXmlElement* lpXmlDir = lpXmlInter->FirstChildElement("Directory") ;
  if ( ! lpXmlDir ) {
    lIntDir = DEFAULT_F77_PROTON_DIR ;
  } else { 
    TiXmlNode* lpXmlDirName = lpXmlDir->FirstChild() ;
    if ( (! lpXmlDirName) || (! lpXmlDirName->ToText()) ) 
      throw TXmlErr("Incorrect interaction directory name");
    lIntDir = lpXmlDirName->Value() ;
  }

  // Maximum time step for interactions
  TiXmlElement* lpXmlStep = lpXmlInter->FirstChildElement("MaxStep_Mpc") ;
  const char* lStep = lpXmlStep->Attribute("value",&_fInteractionTimeStep) ;
  if ( !lStep ) throw TXmlErr("Error getting value of interaction step.") ;
  _fInteractionTimeStep *= Mpc * inv_c_light;

  // Reading interaction tables
  char cha[256] ;
  double a,b ;
  double Lff[81][81],Lpn[81][81],Lnp[81][81] ;
  double deltaE[81] ;
  double lp_int[81],ln_int[81] ;
  double Lmaxp[81],Lmaxn[81] ;
  int nmaxp[81],nmaxn[81] ;
  double Eminpp = 0 ;
  double Eminpn = 0 ;
  double Eminnn = 0 ;
  double Eminnp = 0 ;

  // In this framework, no easy possibility to use the extended pair production table
  string lIntFile = lIntDir+"totppp.out" ;
  fstream totppp(lIntFile.c_str(), ios::in) ;
  if (! totppp.is_open()) {
    throw TCrpErr("Error opening interaction file : " + lIntFile);
  }
  totppp.getline(cha,100) ;
  totppp.getline(cha,100) ;
  for (int i=0; i<=80; i++) {
    totppp >> _fE_part[i] >> _fL_pair[i] ;
    _fE_part[i] /= 1e18 ;
    _fL_pair[i] *= 1e6 ;
  }
  if (!totppp.good()) throw TCrpErr("Error reading totppp.out") ;
  totppp.close() ;
  _fdEtabbin = log10( _fE_part[80]/_fE_part[0] )/80 ;

  lIntFile = lIntDir+"difpp.out" ;
  fstream difpp(lIntFile.c_str(), ios::in) ;
  if (! difpp.is_open()) {
    throw TCrpErr("Error opening interaction file : " + lIntFile);
  }
  difpp.getline(cha,100) ;
  difpp.getline(cha,100) ;
  difpp.getline(cha,100) ;
  for (int i=0; i<=80; i++) {
    for (int j=0; j<=80; j++) difpp >> a >> b >> Lff[i][j] >> Lpn[i][j] >> Lnp[i][j] ;
  }
  if (!difpp.good()) throw TCrpErr("Error reading difpp.out") ;
  difpp.close() ;

  lIntFile = lIntDir+"totpp.out" ;
  fstream totpp(lIntFile.c_str(), ios::in) ;
  if (! totpp.is_open()) {
    throw TCrpErr("Error opening interaction file : " + lIntFile);
  }
  totpp.getline(cha,100) ;
  totpp.getline(cha,100) ;
  totpp.getline(cha,100) ;
  for (int i=0; i<=80; i++) {
    totpp >> a >> _flossp_tot[i] >> _flossn_tot[i] ;
    _flossp_tot[i] *= 1e6 ;
    _flossn_tot[i] *= 1e6 ;
  }
  if (!totpp.good()) throw TCrpErr("Error reading totpp.out") ;
  totpp.close() ;

  for (int i=0; i<=80; i++) {
    for (int j=0; j<=80; j++) {
      if ( j>0 && j<80 ) { // ameliorations simples a faire..
	deltaE[j] = (_fE_part[j+1] - _fE_part[j-1])/2 ;
	Lff[i][j] *= 2 / (_fE_part[j+1]-_fE_part[j-1]) ;
	Lpn[i][j] *= 2 / (_fE_part[j+1]-_fE_part[j-1]) ;
	Lnp[i][j] *= 2 / (_fE_part[j+1]-_fE_part[j-1]) ;
      } else if ( j == 0 ) {
	deltaE[j] = _fE_part[1]-_fE_part[0] ;
	Lff[i][j] /= (_fE_part[1]-_fE_part[0]) ;
	Lpn[i][j] /= (_fE_part[1]-_fE_part[0]) ;
	Lnp[i][j] /= (_fE_part[1]-_fE_part[0]) ;
      } else if ( j == 80 ) {
	deltaE[j] = _fE_part[80]-_fE_part[79] ;
	Lff[i][j] /= (_fE_part[80]-_fE_part[79]) ;
	Lpn[i][j] /= (_fE_part[80]-_fE_part[79]) ;
	Lnp[i][j] /= (_fE_part[80]-_fE_part[79]) ;
      }
    }
  }

  for (int i=0; i<=80; i++) {
    lp_int[i] = 0 ; // improve (put upper)
    ln_int[i] = 0 ;
    for (int j=0; j<=80; j++) {
      lp_int[i] += (Lff[i][j]+Lpn[i][j])*deltaE[j] ;
      ln_int[i] += (Lff[i][j]+Lnp[i][j])*deltaE[j] ;
    }
    _finc_p[i] = 1e6*lp_int[i]/_flossp_tot[i] ;
    _finc_n[i] = 1e6*ln_int[i]/_flossn_tot[i] ;
  }

  for (int i=0; i<=80; i++) {
    Lmaxp[i] = 1e-30 ;
    Lmaxn[i] = 1e-30 ;
    for (int j=0; j<=80; j++) {
      _flossp[i][j] = Lff[i][j]+Lpn[i][j] ;
      _flossn[i][j] = Lff[i][j]+Lnp[i][j] ;
      if ( Lmaxp[i] < _flossp[i][j] ) {
	Lmaxp[i] = _flossp[i][j] ;
	nmaxp[i] = j ;
      }
      if ( Lmaxn[i] < _flossn[i][j] ) {
	Lmaxn[i] = _flossn[i][j] ;
	nmaxn[i] = j ;
      }
    }
    for (int j=0; j<=80; j++) {
      _fLpp[i][j] = Lff[i][j]/Lmaxp[i] ;
      Lpn[i][j] /= Lmaxp[i] ;
      _flossp[i][j] /= Lmaxp[i] ;
      _fLnn[i][j] = Lff[i][j]/Lmaxn[i] ;
      Lnp[i][j] /= Lmaxn[i] ;
      _flossn[i][j] /= Lmaxn[i] ;
    }
    for (int j=0; j<=79; j++) {
      if (_fLpp[i][j] < 1e-7 && _fLpp[i][j+1] >= 1e-7) Eminpp = _fE_part[j] ;
      if (Lpn[i][j] < 1e-7 && Lpn[i][j+1] >= 1e-7) Eminpn = _fE_part[j] ;
      if (_fLnn[i][j] < 1e-7 && _fLnn[i][j+1] >= 1e-7) Eminnp = _fE_part[j] ;
      if (Lnp[i][j] < 1e-7 && Lnp[i][j+1] >= 1e-7) Eminnn = _fE_part[j] ;
    }
    _fEminp[i] = min(Eminpp, Eminpn) ; 
    if ( _fEminp[i] == 0 ) _fEminp[i] = _fE_part[0] ;
    _fEminn[i] = min(Eminnp, Eminnn) ; 
    if ( _fEminn[i] == 0 ) _fEminn[i] = _fE_part[0] ;
  }

  _finc_p[0] = 0 ;
  _finc_n[0] = 0 ;

  // Neutron decay table : energy distribution of outcoming electron
  double lQ = neutron_mass_c2 - proton_mass_c2 ;
  int lSize = 100 ;
  double* lEnergies = new double[lSize] ;
  double* lDensity = new double[lSize] ;
  for (int i=0; i<lSize; i++) {
    lEnergies[i] = electron_mass_c2 + i*(lQ-electron_mass_c2)/(lSize-1) ;
    lDensity[i] = lEnergies[i]*
      sqrt(max(lEnergies[i]*lEnergies[i]-electron_mass_c2*electron_mass_c2,0.))*
      (lQ-lEnergies[i])*(lQ-lEnergies[i]) ;
  }
  _fpRandDistriNeutronDecay = new RandGeneral(lDensity,lSize,0) ;
  delete[] lEnergies ;
  delete[] lDensity ;

  // Switches interactions on/off
  TiXmlElement* lpXmlCheckRedshift = lpXmlInter->FirstChildElement("NoRedshift") ;
  if (!lpXmlCheckRedshift) {
    _fRedshiftFlag = 1 ;
  } else _fRedshiftFlag = 0 ;
  TiXmlElement* lpXmlCheckPhoton = lpXmlInter->FirstChildElement("SecondaryPhotons") ;
  if ( lpXmlCheckPhoton ) 
    throw TCrpErr("No secondary photons with Proton-F77.") ;
  TiXmlElement* lpXmlCheckNu = lpXmlInter->FirstChildElement("SecondaryNeutrinos") ;
  if ( lpXmlCheckNu )
    throw TCrpErr("No secondary neutrinos with Proton-F77.") ;
  TiXmlElement* lpXmlCheckPairP = lpXmlInter->FirstChildElement("NoPairProd") ;
  if ( !lpXmlCheckPairP ) {
    _fPairProdFlag = 1 ;
  } else _fPairProdFlag = 0 ;
  TiXmlElement* lpXmlCheckPionP = lpXmlInter->FirstChildElement("NoPionProd") ;
  if ( !lpXmlCheckPionP ) {
    _fPionProdFlag = 1 ;
  } else _fPionProdFlag = 0 ;

  if (!_fPionProdFlag && !_fPairProdFlag) 
    throw TXmlErr("Proton-F77 but no interactions are left on.") ;

}

TBasicPInteractions::~TBasicPInteractions() {

  delete _fpRandDistriNeutronDecay ;

}
