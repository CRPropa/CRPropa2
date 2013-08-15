/**
   @file    clustergas.h
   @author  Luca Maccione, luca.maccione@desy.de
   @brief   Class describing a typical gas distribution in clusters, as used in Kotera et al. Astrophys.J. 707 (2009) 370-386
*/

#ifndef _CLUSTERGAS_H_
#define _CLUSTERGAS_H_

#include "gas.h"
#include "xmlparam.h"
#include "xmlextract.h"
#include "crp_err.h"

#include <vector>
#include <string>
#include <fstream>
#include <math.h>

/**
   @class TClusterGas
   @brief Gas density

*/

using namespace std;

class TClusterGas: public TGas, public TXmlParam {

 public:
  TClusterGas(const char*);
  ~TClusterGas() { }

  double Rmax() const { return _fRarray.back(); }
  int    Nr()   const { return _fRarray.size(); }
  // vector<double> Position() const { return _fRarray; }
  double Position(int i) const {return _fRarray.at(i); }
  // vector<double> Density() const { return _fDensity; }
  // double Density(int i) const { return _fDensity.at(i); }
  virtual double Density(double);

 protected:
  vector<double>_fRarray;
  vector<double>_fDensity;

};

#endif
