/**
   @file    variableinfrared.h
   @author  Luca Maccione, luca.maccione@desy.de
   @brief   Class describing a variable infrared background, with simple redshift evolution
*/


#ifndef _VARIABLEINFRARED_H_
#define _VARIABLEINFRARED_H_

#include "xmlparam.h"
#include "xmlextract.h"
#include "crp_err.h"
#include "IR.h"
#include "PhotonBackground.h"
#include "photonspectrum.h"

#include <vector>
#include <string>
#include <fstream>
#include <math.h>
#include <algorithm>

//#define DEFAULT_ZMAX_IR 5.
/**< Default value for the redshift where the IR is generated. It can be redefined in the configuration file. */

/**
   @class TVariableInfrared
   @brief Infrared background
*/

using namespace std;

class TVariableInfrared : public IR, public TXmlParam {
 public:
  TVariableInfrared(const char*, bool) ;
  virtual ~TVariableInfrared() { } 
    
  virtual double Zmax() { return _fZmax; }
  virtual double GetPhotonDensity(double x, 
				  double y, 
				  double z, 
				  double redshift,
				  double PhotonEnergy /* GeV */);

  /**< Maximum redshift of the IR background. */
  virtual double Epsmin() { return _fEnergy.front(); }
  virtual double Epsmax() { return _fEnergy.back(); }
  virtual TPhotonSpectrum Spectrum(double x,
				    double y = -1,
				    double z = -1);
  virtual vector<double> GetE() { return _fEnergy; }
  virtual vector<double> GetX() { return _fPositionX; }
  virtual vector<double> GetY() { return _fPositionY; }
  virtual vector<double> GetZ() { return _fPositionZ; }
  virtual vector<double> GetR() { return this->GetX(); }

  virtual int NE() { return _fNE; }
  virtual int Nx() { return _fNx; }
  virtual int Ny() { return _fNy; }
  virtual int Nz() { return _fNz; }
  virtual int Nr() { return this->Nx(); }

 private:
  double _fZmax ;
  int _fNx;
  int _fNy;
  int _fNz;
  int _fNE;
  int _fNdim;
  vector<double> _fPositionX;
  vector<double> _fPositionY;
  vector<double> _fPositionZ;
  vector<double> _fEnergy;
  vector< vector<double> > _fEnergyDensity;
  bool oned;
};


void readstring(const string&, vector<double>&);
void readstring(const string&, vector<double>&, vector< vector<double> >&);

#endif
