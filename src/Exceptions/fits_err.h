
/**
	@file		fits_err.h
	@author		Eric Armengaud, armengau@in2p3.fr
	@brief		Exceptions taking place in cfitsio interface
*/

#ifndef _FITS_ERR_H_
#define _FITS_ERR_H_

#include "crp_err.h"

/**
   @class TFitsErr
   @brief Class inherited from TCrpErr to handle errors occuring from the CFITSIO module

   This simply allows to let the user know that the error is due to CFITSIO, and also to know the CFITSIO error code associated. A list of these error codes can be found for example on :
   http://heasarc.gsfc.nasa.gov/docs/software/fitsio/quick/node26.html
*/

class TFitsErr : public TCrpErr {
 public:
  TFitsErr() { }
  TFitsErr(const char *aMessage) { 
    std::cerr << "CFITSIO error." << std::endl ;
    std::cerr << aMessage << std::endl;
    exit(ERR_FITS);
  } /**< Prints a message to error stream, and leave the simulation */
  TFitsErr(std::string aMessage) {
    std::cerr << "CFITSIO error. " << std::endl ;
    std::cerr << aMessage << std::endl;
    exit(ERR_FITS);
  } /**< Prints a message to error stream, and leave the simulation */
  TFitsErr(int status) {
    std::cerr << "Fits output error. " << std::endl ;
    std::cerr << "CFITSIO error code is: " << status << std::endl ;
    exit(ERR_FITS) ;
  } /**< Prints the CFITSIO error code to error stream, and leave the simulation */

  ~TFitsErr() { }
} ;

#endif
