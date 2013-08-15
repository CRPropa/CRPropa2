/**
   @file   xmlparam.h
   @author Tristan Beau, beau@in2p3.fr
   @brief  Generic class to get parameters from a xml file
*/

#ifndef _XMLPARAM_H_
#define _XMLPARAM_H_

#include "xmlextract.h"

/**
   @class TXmlParam
   @brief Generic class that contains a TXmlExtract object. Every class whose objects need to access 
the XML configuration file (for instanciation for example) should inherit from this class.
*/

class TXmlParam {
 protected:
  TXmlParam() {_fpXmlExtract = NULL;}
  TXmlParam(const char* aFileName) { _fpXmlExtract = new TXmlExtract(aFileName); }
  ~TXmlParam() { delete _fpXmlExtract; }

  TXmlExtract &XmlExtract() const { return *_fpXmlExtract; } 
  /**< Returns a pointer to the XmlExtract member.*/

 private:
  TXmlExtract *_fpXmlExtract;

};

#endif
