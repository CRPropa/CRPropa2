/**
   @class IR

   @author  Luca Maccione (maccioneluca@googlemail.com) \n 
            Nils Nierstenhoefer (nierstenhoefer@physik.uni-wuppertal.de)


   @brief   Returns photon density (h_bar=c=1) in units of GeV^2   

   A "cmb-like" scaling of IR photon density with redshift is used (s. CPRropa paper equation 12)
   
   \f$
   
   n[\epsilon, z] =  
   \left\{
   \begin{array}{l l l}
   
   & (1+z)^2 \cdot n[\frac{\epsilon}{1+z}, z=0] & z \leq z_{\rm{b}}
   \\
   & 0 & \rm{otherwise.} 
   
   \end{array}
   \right\}  
   \f$

   Here, \f$ z_{\rm{b}}\sim 1\f$ is used.
*/

#ifndef _IR_H_
#define _IR_H_
//#define DEFAULT_ZMAX_IR 5

#include "typedclass.h"
#include "PhotonBackground.h"

class IR : public PhotonBackground {
 public:
  
 IR():_fZmax(DEFAULT_ZMAX_IR){ 
    fEpsilonZero=4.6e-12; 
  }
  
  IR(double _fZmax_) : _fZmax(_fZmax_){ 
    fEpsilonZero=4.6e-12;
  }
    
    virtual ~IR() {
      //std::cout<<"fEpsilonZero set to "<<fEpsilonZero<<std::endl;
    }
    
  virtual double GetPhotonDensity(double x, 
				  double y, 
				  double z, 
				  double redshift,
				  double PhotonEnergy);
  virtual double Zmax() { return _fZmax; };
  virtual int    NE() { return 0;}
  virtual int Nx() { return 0;}
  virtual int Ny() { return 0;}
  virtual int Nz() { return 0;}
  virtual int Nr() { return 0;}
  virtual double Epsmin() { return 0;}
  virtual double Epsmax() { return 0;}
  virtual std::vector<double> GetE() { return vector<double>();}
  virtual std::vector<double> GetX() { return vector<double>();}
  virtual std::vector<double> GetY() { return vector<double>();}
  virtual std::vector<double> GetZ() { return vector<double>();}
  virtual std::vector<double> GetR() { return vector<double>();}
  virtual TPhotonSpectrum Spectrum(double x, double y = -1, double z = -1) { return TPhotonSpectrum();}

 protected :
  double _fZmax;
};





#endif

