/**
   @class PhotonBackground

   @author  Luca Maccione (maccioneluca@googlemail.com),\n
            Nils Nierstenhoefer (nierstenhoefer@physik.uni-wuppertal.de)

   @brief Prototype class to handle photon fields. 

   <strong> The Case of "Constant" Photon Fields </strong>

   Although many CRPropa routines are in principle able to handle variable photon fields
   \f$n(z, \vec{x})\f$, CRPropa v1.5 is restricted to constant photon fields. This is mainly
   due to corresponding limitations in the actual propagation algorithm. Thus, the mean free path
   is always calculated using the photon density at an redshift z=0. The redshift dependence then
   comes in by assuming "CMB-like" evolution which gives rise to the following scaling relation
   for the reaction rates

   \f$ \lambda^{-1}=(1+z)^3 \cdot \lambda[(1+z)E,z=0]^{-1}.\f$.

   <strong> Todo </strong>

   As some functionality of the proton propagation of CRPropa 1.3 have been generalized to the case
   of nuclei some parts of CRPropa 1.5 (e.g. pion production) still rely on the photon fields defined
   in the fortran code given in <CODE>$CRPropaSrc/sophia/src/sophia_interface.f</CODE>. I.e. a change
   of the C++ code documented here will have no effect in this case! This is a historical grown difficulty 
   that should be sorted out in a future release.
   
   The Get GetEpsilonZero() member was only introduced for the calculation of the reaction rate 
   for nuclei photo disintegration - integral substitution \f$t=ln(\epsilon/\epsilon_{0})\f$. Thus, it should 
   probably be moved to the TALYSMeanFreePath class.
*/

#ifndef _PhotonBackground_H_
#define _PhotonBackground_H_

#include "typedclass.h"
#include "photonspectrum.h"
#include <vector>
#include <TGraph.h>

#define DEFAULT_ZMAX_IR 5

class PhotonBackground : public TTypedClass {
 public:
  virtual double GetPhotonDensity(double x,
				  double y,
				  double z,
				  double redshift,
				  double PhotonEnergy)=0;
  /**< Returns the photon density in units of \f$\rm{GeV}^2\f$ */
  virtual double GetEpsilonZero() {return fEpsilonZero;} 
  /**< For the photo disintegration a  integral substitution \f$t=ln(\epsilon / \epsilon_{0})\f$ was performed. This routine provides the corresponding \f$ \epsilon_{0} \f$ value. */ 
  virtual double Zmax() { return 0;}
  virtual int    NE()  { return 0;}
  virtual int Nx() { return 0;}
  virtual int Ny() { return 0;}
  virtual int Nz() { return 0;}
  virtual int Nr() { return 0;}
  virtual double Epsmin() { return 0;}
  virtual double Epsmax() { return 0;}
  virtual std::vector<double> GetE() { return vector<double>();}
  virtual std::vector<double> GetR() { return vector<double>();}
  virtual std::vector<double> GetX() { return vector<double>();}
  virtual std::vector<double> GetY() { return vector<double>();}
  virtual std::vector<double> GetZ() { return vector<double>();}
  virtual TPhotonSpectrum Spectrum(double x, double y = -1, double z = -1) { return TPhotonSpectrum();}
  
  virtual ~PhotonBackground() { }
  PhotonBackground() //:
    //   fEpsilonZero(0)
    {}

  double fEpsilonZero;

 protected:
     
 private:
 


};


//Plot function (not a class member).
TGraph* RootPlotSpectrum(PhotonBackground& BGrd,
			 double Emin,
			 double Emax,
			 int EBins,
			 double x,
			 double y,
			 double z,
			 double redshift);



#endif

