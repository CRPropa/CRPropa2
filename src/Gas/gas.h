
/**
   @file    gas.h
   @author  Luca Maccione, luca.maccione@desy.de
   @brief   Virtual class describing IG gas
*/

#ifndef _GAS_H_
#define _GAS_H_

#include "typedclass.h"
#include "vector3d.h"
#include <fstream>
#include "units.h"
#include "crp_err.h"
#include "xml_err.h"
#include "fits_err.h"
#include <vector>

using namespace std;

/**
   @class TGas
   @brief Virtual class for IG gas. Currently there is only one implementation.

*/

class TGas : public TTypedClass {
 public:
	
  virtual ~TGas() { }
  virtual double Rmax() const { return 0.; }
  virtual int Nr() const { return 0; }
  virtual int Nx() const { return 0; }
  virtual int Ny() const { return 0; }
  virtual int Nz() const { return 0; }
  virtual double Stepsize() const { return 0;}
  //  virtual vector<double> Position() const { vector<double> result; return result;}
  virtual double Position(int) const { return 0.0;}
  // virtual vector<double> Density() { vector<double> result; return result; }
  //virtual double Density(int) const { return 0.0; }
  virtual double Density(TVector3D) { return 0.0; }
  virtual double Density(double) { return 0.0; }
  virtual double Gasmax() const { return 0; }
  virtual TVector3D Origin() const { return TVector3D(); }
};

#endif
