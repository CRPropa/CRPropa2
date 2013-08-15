/**
   @class   CMB

   @author  Nils Nierstenhoefer (nierstenhoefer@physik.uni-wuppertal.de)

   @brief   Returns photon density (h_bar=c=1) in units of GeV^2   

   The photon density scales with redshift according to (s CPRropa paper equation 12)
   
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

   Here, \f$ z_{\rm{b}}\sim 1100\f$ is the redshift of decoupling.
*/

#ifndef _CMB_H_
#define _CMB_H_

#include "typedclass.h"
#include "PhotonBackground.h"

class PhotonBackground;

class CMB : public PhotonBackground {
 public:
  CMB();
  ~CMB();
  virtual double GetPhotonDensity(double x, 
				  double y, 
				  double z, 
				  double redshift,
				  double PhotonEnergy);

  virtual double GetEpsilonZero() {return fEpsilonZero;}

 protected :
  //double fEpsilonZero;

 private :

};





#endif

