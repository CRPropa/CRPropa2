/**
   @file    magfield.h
   @author  Eric Armengaud, armengau@in2p3.fr
   @brief   Virtual class describing any magnetic field in the simulation
*/

#ifndef _MAGFIELD_H_
#define _MAGFIELD_H_

#include "typedclass.h"
#include "vector3d.h"
#include <fstream>
#include "units.h"
#include "crp_err.h"
#include "xml_err.h"
#include "fits_err.h"

using namespace std;

/**
   @class TMagField
   @brief Virtual class for the magnetic fields.

   There are currently four implementations of this class: the TNullMagField and TUniformMagField are trivial. The TGridField allows to implement any field that is defined on a regular 3-dimensionnal grid. The TField1D describes the magnetic field perpendicular to the line of sight in 1-D simulations, and is useful to compute the synchrotron radiation in secondary electromagnetic cascades.

*/

class TMagField : public TTypedClass {
 public:
	
  virtual ~TMagField() {}

  virtual TVector3D* B() const { return NULL; }
  virtual int Nx() const { return 0; }
  virtual int Ny() const { return 0; }
  virtual int Nz() const { return 0; }
  virtual double Stepsize() const { return 0; }
  virtual TVector3D Origin() const { return TVector3D(); }

  virtual TVector3D getField(TVector3D) const { return TVector3D(); }

  virtual double getPositions(int) const { return 0; }
  virtual double getFieldValues(int) const { return 0; }
  virtual double Xmax() const { return 0; }

  virtual double Bmax() const { return 0; }
  virtual int writeMagField(string) {return 0;};

};

#endif
