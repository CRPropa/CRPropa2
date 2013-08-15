/**
   @file    kolmogoroffmagfield.h
   @author  Eric Armengaud, armengau@in2p3.fr
   @brief   Class describing a magnetic field with Kolmogoroff spectrum
*/

#ifndef _KOLMOGOROFFMAGFIELD_H_
#define _KOLMOGOROFFMAGFIELD_H_

#include "gridfield.h"
#include "xmlparam.h"
#include "crp_err.h"

#include <math.h>
#include <string>

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGauss.h"
#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

/**
   @class TKolmogoroffMagField
   @brief Statistically homogeneous, turbulent magnetic field implemented on a regular grid.
   
   Before the simulation, a field configuration is computed with random amplitudes and phases, on a regular grid. This field therefore inherits from the TGridField class. It is defined in Fourier space by a random gaussian amplitude B(k) with a single power law distribution in a given spectral range.

*/

class TKolmogoroffMagField : public TGridField, public TXmlParam {
 public:
  TKolmogoroffMagField(const char*) ;
  ~TKolmogoroffMagField() ;

  void Output(int, int, int);

  double SpectralIndex() const { return _fSpectralIndex; }
  /**< Spectral index of the power law distribution for the amplitude B^2(k). The case n=-11/3 corresponds to a Kolmogoroff spectrum, but in fact any other spectral index can be specified. */
  double Kmin() const { return _fKmin; }
  /**< Minimum k value for the power law distribution of B(k). */
  double Kmax() const { return _fKmax; }
  /**< Maximum k value for the power law distribution of B(k). */
  double Bmax() const { return _fBmax; }
  /**< Maximum value of the field on the grid. */

  void Fourn(float[], unsigned long[], int, int) ;
  /**< Routine used to compute the field in real space. Inspired from Num. Recipes. */

 protected:
  double _fSpectralIndex ;
  double _fKmin, _fKmax ;
  double _fBmax ;
  string _fOutputFileName;
  ofstream _fMagFieldStream ;

};

#endif
