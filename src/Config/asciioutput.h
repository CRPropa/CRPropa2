
/**
   @file   asciioutput.h
   @author Eric Armengaud, armengau@in2p3.fr
   @brief  Class describing ascii format output
*/

#ifndef _ASCIIOUTPUT_H_
#define _ASCIIOUTPUT_H_

#include "outputdata.h"
#include "xmlparam.h"
#include "units.h"
#include "vector3d.h"

#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

/**
   @class TAsciiOutput
   @brief Class allowing to write the output into an ASCII file

   The output file is open and a header is written when this class is instanciated.
   It is then filled during propagation with the functions Add1/3DEvent/Traj/Shower/Neutrino,
   and closed when the object is destructed.
   The arguments of Add(..) functions must have CRPropa units.
   For secondary neutrinos and/or electromagnetic cascades, two or three files are open and concatenated in the end.
   In the case of electromagnetic cascades, a first line is written with the list of energies corresponding 
   to the output spectrum (in EeV).

   The type is OUTPUT_ASCII.

*/

class TAsciiOutput : public TOutputData, public TXmlParam {
 public:
  TAsciiOutput(const char*, string, bool aForceOverwrite = false) ; 
  /**< Opens the files. aForceOverwrite allows to overwrite the output even if it already exists. */
  ~TAsciiOutput() ;
  
  void Add1DEvent(int, double, double, double, double, double, int) ; 
  /**< Arguments are : type - initial position - initial energy - time - energy - initial Type */
  void Add3DEvent(int, TVector3D, TVector3D, double, TVector3D, TVector3D, int) ;
  /**< Arguments are : type - initial position - initial momentum - time - position - momentum - initial Type */

  void Add1DTraj(int, double, double, double, int) ;
  /**< Arguments are : type - time - position - energy - Initial Type */
  void Add3DTraj(int, double, TVector3D, TVector3D, double, int) ;
  /**< Arguments are : type - time - position - momentum - energy - Initial Type */
  void Add3DTraj(int, double, TVector3D, TVector3D, double,
		 double, double, int) ;
  /**< Arguments are : type - time - position - energy - momentum - local B field - local source density - Initial Type This function
   is not used currently. */

  void Add1DShower(string, double, double, double, double, vector<double>) ;
  /**< Arguments are : origin - source position - initial position - energy - source energy - spectrum */
  void Add3DShower(string, TVector3D, TVector3D, TVector3D, TVector3D, double,
		   vector<double>, int) ;
  /**< Arguments are : origin - source position - initial position - position - momentum - source energy - spectrum - Initial Type */

  void Add1DNeutrino(int, double, double, double, double, double, int ) ;
  /**< Arguments are : type - source position - initial position - initial energy - energy - source energy - initial Type */
  void Add3DNeutrino(int, TVector3D, TVector3D, TVector3D, TVector3D, double, int) ;
  /**< Arguments are : type - source position - initial position - position - momentum - source energy - initial Type */

 private:
  ofstream _fDataStream ;
  ofstream _fShowerStream ;
  ofstream _fNeutrinoStream ;
};

#endif
