/**
   @file   outputdata.h
   @author Eric Armengaud, armengau@in2p3.fr
   @brief  Class describing the interface for output data
*/

#ifndef _OUTPUTDATA_H_
#define _OUTPUTDATA_H_

#include "typedclass.h"
#include "xmlparam.h"
#include <vector>
#include <fstream>
#include "units.h"
#include "crp_err.h"
#include "xml_err.h"

#include "fitsio.h"

#include "vector3d.h"

using namespace std;

/**
   @class TOutputData
   @brief Virtual class describing the interface of output data.
*/

class TOutputData : public TTypedClass {
 public:
  virtual ~TOutputData() {}

  virtual void Add1DEvent(int, double, double, double, double, double, int) {}
  virtual void Add3DEvent(int, TVector3D, TVector3D, double, TVector3D, TVector3D, int) {}
  virtual void Add1DTraj(int, double, double, double, int) {}
  virtual void Add3DTraj(int, double, TVector3D, TVector3D, double, int) {}
  virtual void Add3DTraj(int, double, TVector3D, TVector3D,
			 double, double, double, int) {}
  virtual void Add1DShower(string, double, double, double, double,
			     vector<double>) {}
  virtual void Add3DShower(string, TVector3D, TVector3D, TVector3D, TVector3D, double,
			   vector<double>, int) {}
  virtual void Add1DNeutrino(int, double, double, double, double, double, int) {}
  virtual void Add3DNeutrino(int, TVector3D, TVector3D, TVector3D, TVector3D, double, int) {}

 protected:
  string _fFileName ;
};

#endif
