
/**
   @file    photoninteractions.cc
   @author  Eric Armengaud, armengau@in2p3.fr
   @brief   Implementation of the TPhotonInteractions class. See the .h file
*/

#include "photoninteractions.h"

TPhotonInteractions::TPhotonInteractions(const char* aFileName) : TXmlParam(aFileName) {

  SetType(INTERACTION_PHOTON);
  TiXmlElement* lpXmlInter = XmlExtract().GetElement("Interactions") ;

  // Directory for shower tables
  TiXmlElement* lpXmlShowerDir = lpXmlInter->FirstChildElement("ShowerTableDirectory") ;
  if ( ! lpXmlShowerDir ) {
    //      throw TXmlErr("Incorrect photon shower table directory");
    _fShowerTableDir = DEFAULT_SHOWER_DIR ;
  } else {
    TiXmlNode* lpXmlShowerDirName = lpXmlShowerDir->FirstChild() ;
    if ( (! lpXmlShowerDirName) || (! lpXmlShowerDirName->ToText()) )
      throw TXmlErr("Incorrect photon shower table directory name");
    _fShowerTableDir = lpXmlShowerDirName->Value() ;
  }

  // Cascade cut
  TiXmlElement* lpXmlCheckCutcascade = lpXmlInter->FirstChildElement("CutCascade_MagField") ;
  if ( !lpXmlCheckCutcascade ) {
    _fCutcascadeFlag = 0 ;
  } else {
    const char* lDum = lpXmlCheckCutcascade->Attribute("value",&_fCutcascadeFlag);
    if ( lDum == NULL ) _fCutcascadeFlag = 1 ;
  }

  // Pair production table : in the same directory as fShowerTableDir.
  string lPairProdTable = _fShowerTableDir+"/pair_spectrum_cmbir.table";
  fstream pptable(lPairProdTable.c_str(),ios::in);
  if (! pptable.is_open()) throw TCrpErr("Error opening pair production table : "+lPairProdTable);
  char cha[256]; double a;
  pptable.getline(cha,200);
  pptable.getline(cha,200);
  pptable.getline(cha,200);
  for (int i=0; i<NBINS_PAIRPROD; i++) {
    for (int j=0; j<NUM_MAIN_BINS; j++) {
      pptable >> a;
      _fPairProdTable[i][j]=a;
    }
  }
  if (!pptable.good()) throw TCrpErr("Error reading pairprod_table.txt") ;
  pptable.close() ;
  double lEpmin=pow(10.,EPMIN_PAIRPROD);
  double lEpmax=pow(10.,EPMAX_PAIRPROD);
  for (int i=0; i<NBINS_PAIRPROD; i++)
    _fNucleonEnergy_PP[i]=lEpmin*pow(lEpmax/lEpmin,double(i)/double(NBINS_PAIRPROD));
  
}

TPhotonInteractions::~TPhotonInteractions() {

}

double TPhotonInteractions::PairProdSpec(double aEp_eV, int aBin_Ee) const {
  double lInterpolValue=0;
  if (aEp_eV < _fNucleonEnergy_PP[0]) lInterpolValue=0;
  else if (aEp_eV > _fNucleonEnergy_PP[NBINS_PAIRPROD-1]) lInterpolValue=0;
  else {
    int lBin_Ep=0;
    while (aEp_eV > _fNucleonEnergy_PP[lBin_Ep]) {
      lBin_Ep++;
      if (lBin_Ep == NBINS_PAIRPROD-2) break;
    }
    double alpha=(aEp_eV-_fNucleonEnergy_PP[lBin_Ep])/(_fNucleonEnergy_PP[lBin_Ep+1]-_fNucleonEnergy_PP[lBin_Ep]);
    lInterpolValue = (1-alpha)*_fPairProdTable[lBin_Ep][aBin_Ee] + alpha*_fPairProdTable[lBin_Ep+1][aBin_Ee];
  }
  
  return lInterpolValue;
}
