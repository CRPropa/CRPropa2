
/**
	@file		xml_err.h
	@author		Eric Armengaud, armengau@in2p3.fr
	@brief		Exceptions taking place while parsing XML files
*/

#ifndef _XML_ERR_H_
#define _XML_ERR_H_

#include "crp_err.h"

/**
   @class TXmlErr
   @brief Class inherited from TCrpErr to handle errors occuring from the TinyXML module

   This simply allows to let the user know that the error is due to TinyXml
*/

class TXmlErr : public TCrpErr {
 public:
  TXmlErr() { }
  TXmlErr(const char *aMessage, bool isFatal=true) {
    if ( isFatal )  { 
      std::cerr << "XML file parsing error. " << std::endl ;
      std::cerr << aMessage << std::endl;
      exit(ERR_XML);
    }
  }  /**< Prints a message to error stream, and leave the simulation if the error flag isFatal is true */
  TXmlErr(std::string aMessage, bool isFatal=true) {
    if ( isFatal ) {
      std::cerr << "XML file parsing error. " << std::endl ;
      std::cerr << aMessage << std::endl;
      exit(ERR_XML);
    }
  }  /**< Prints a message to error stream, and leave the simulation if the error flag isFatal is true */
  ~TXmlErr() { }
} ;

#endif
