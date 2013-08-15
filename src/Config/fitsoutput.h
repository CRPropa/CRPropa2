
/**
   @file   fitsoutput.h
   @author Eric Armengaud, armengau@in2p3.fr
   @brief  Class describing Fits format output
*/

#ifndef _FITSOUTPUT_H_
#define _FITSOUTPUT_H_

#include "fitsio.h"

#include "fits_err.h"
#include "outputdata.h"
#include "xmlparam.h"

#include <fstream>
#include "units.h"

using namespace std;

/**
   @class TFitsOutput
   @brief Class allowing to write the output into a FITS file

   The output file is open and a header is written when this class is instanciated.
   It is then filled during propagation with the functions Add1/3DEvent/Traj/Shower/Neutrino,
   and closed when the object is destructed.
   The arguments of Add(..) functions must have CRPropa units.
   If secondaries are recorded, several HDUs are filled.
   In the case of electromagnetic cascades, a first line is written with the list of energies corresponding 
   to the output spectrum (in EeV).

   The type is OUTPUT_FITS.

*/

class TFitsOutput : public TOutputData, public TXmlParam {
 public:
  TFitsOutput(const char*, string, bool aForceOverwrite = false) ;
  /**< Opens the files. aForceOverwrite allows to overwrite the output even if it already exists. */
  ~TFitsOutput() ;

  void Add1DEvent(int, double, double, double, double, double, int) ;
  /**< Arguments are : type - initial position - initial energy - time - energy - Initial Type*/
  void Add3DEvent(int, TVector3D, TVector3D, double, TVector3D, TVector3D, int) ;
  /**< Arguments are : type - initial position - initial momentum - time - position - momentum - Initial Type */

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
  void Add3DShower(string, TVector3D, TVector3D, TVector3D, TVector3D,
		   double, vector<double>, int) ;
  /**< Arguments are : origin - source position - initial position - position - momentum - source energy - spectrum - Initial Type */

  void Add1DNeutrino(int, double, double, double, double, double, int) ;
  /**< Arguments are : type - source position - initial position - initial energy - energy - source energy - Initial Type */
  void Add3DNeutrino(int, TVector3D, TVector3D, TVector3D, TVector3D, double, int) ; /**< Arguments are : type - source position - initial position - position - momentum - source energy - Initial Type */

 private:
  fitsfile *_fpFitsStream ;
  int _fStatus ; /**< Status of the _fpFitsStream */
  int _fNbFields ; /**< Number of fields in the table of _fpFitsStream */
  long _fCurrentRow ; /**< Current row for the _fpFitsStream */

  fitsfile *_fpShowerFitsStream ;
  int _fShowerStatus ; /**< Status of the _fpShowerFitsStream */
  int _fShowerNbFields ; /**< Number of fields in the table of _fpShowerFitsStream */
  long _fShowerCurrentRow ; /**< Current row for the _fpShowerFitsStream */

  fitsfile *_fpNuFitsStream ;
  int _fNuStatus ; /**< Status of the _fpNuFitsStream */
  int _fNuNbFields ; /**< Number of fields in the table of _fpNuFitsStream */
  long _fNuCurrentRow ; /**< Current row for the _fpNuFitsStream */

};

#endif
