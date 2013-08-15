
/**
   @file   rootoutput.h
   @author Eric Armengaud, armengau@in2p3.fr
   @brief  Class describing Root format output
*/

#ifndef _ROOTOUTPUT_H_
#define _ROOTOUTPUT_H_

// Compile that if you have root
#include "sysdep.h"
//#ifdef HAVE_TROOT_H

#include "outputdata.h"
#include "xmlparam.h"
#include "units.h"
#include "vector3d.h"

#include "TNtuple.h"
#include "TFile.h"
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

class TRootOutput : public TOutputData, public TXmlParam {
 public:
  TRootOutput(const char*, std::string, bool aForceOverwrite = false) ;
  ~TRootOutput() ;

  void Add1DEvent(int, double, double, double, double, double, int) ;
  /**< Arguments are : type - initial position - initial energy - time - energy -Initial_Type */
  void Add3DEvent(int, TVector3D, TVector3D, double, TVector3D, TVector3D, int) ;
  /**< Arguments are : type - initial position - initial momentum - time - position - momentum -Initial_Type */

  void Add1DTraj(int, double, double, double, int) ;
  /**< Arguments are : type - time - position - energy - Initial_Type */
  void Add3DTraj(int, double, TVector3D, TVector3D, double, int) ;
  /**< Arguments are : type - time - position - momentum - energy - Initial Type */
  void Add3DTraj(int, double, TVector3D, TVector3D, double,
		 double, double, int) ;
  /**< Arguments are : type - time - position - energy - momentum - local B field - local source density - Initial Type This function
   is not used currently. */

  void Add1DShower(string, /*int,*/ double, double, double, double, vector<double>) ;
  /**< Arguments are : origin - InitParticleType, source position - initial position - energy - source energy - spectrum */
  void Add3DShower(string, TVector3D, TVector3D, TVector3D, TVector3D, double,
		   vector<double>, int) ;
  /**< Arguments are : origin - source position - initial position - position - momentum - source energy - spectrum -Initial_Type */

  void Add1DNeutrino(int, double, double, double, double, double, int) ;
  /**< Arguments are : type - source position - initial position - initial energy - energy - source energy - Initial_Type */
  void Add3DNeutrino(int, TVector3D, TVector3D, TVector3D, TVector3D, double, int) ;
  /**< Arguments are : type - source position - initial position - position - momentum - source energy - Initial Type */

 private:
  TFile* _fRootFile;
  TNtuple* _fEventNtuple;
  TTree* _fShowerTree;
  TNtuple* _fNeutrinoNtuple;
// Parameters for the EM shower TTree
  int _fParticle_Type ;
  int _fInit_Particle_Type;
  //int _fCascadeOrigin ;
  char _fCascadeOriginString[11];
  double _fSourcePosition_X_Mpc ;
  double _fSourcePosition_Y_Mpc ;
  double _fSourcePosition_Z_Mpc ;
  double _fInitialPosition_X_Mpc ;
  double _fInitialPosition_Y_Mpc ;
  double _fInitialPosition_Z_Mpc ;
  double _fPosition_X_Mpc ;
  double _fPosition_Y_Mpc ;
  double _fPosition_Z_Mpc ;
  double _fEnergy_EeV ;
  double _fMomentum_theta ;
  double _fMomentum_phi ;
  double _fSourceEnergy_EeV ;
  float _fSpectrumEnergy_EeV[170] ;
  float _fSpectrum[170] ;

};

//#endif
#endif

