/**
   @file   xmlextract.h
   @author Tristan Beau, beau@in2p3.fr
   @brief  Definition of the TXmlExtract class which extract scalar or vector from xml file. It uses currently the tinyxml package but could be written using another xml parser. 
*/

#ifndef _XMLEXTRACT_H_
#define _XMLEXTRACT_H_

#include "tinyxml.h"
#include <fstream>
#include "xml_err.h"
#include <stdlib.h>
#include <limits.h>
using namespace std;

enum XmlExtractError { 
  ERR_LOAD_FILE = 0x21
};

string GetAttribute(const TiXmlElement*, const char*); /// Tests and returns an attribute for a TiXmlElement

/**
   @class TXmlExtract
   @brief Class containing a pointer to a TiXmlDocument object, and simple functions to check and extract basic elements

   The TiXmlExtract is currently a simple interface between CRPropa and tinyxml objects.
*/

class TXmlExtract {
 public:
  TXmlExtract(const char *)  ;
  ~TXmlExtract() ;
  
  TiXmlElement* GetElement(const char*) const ; /**< Extracts an XML element of a given name from the document */
  int GetInt(const char*) const ; /**< Extracts the integer I for the following XML tag : <Arg value = I /> */
  unsigned long GetULong(const char*) const ; /** Extracts the unsigned long L for the following XML tag :
						  <Arg value = L /> */
  double GetDouble(const char*) const ; /** Extracts the double D for the following XML tag :
					    <Arg value = D /> */
  string GetString(const char*) const ; /** Extracts the string S for the following XML tag :
					    <Arg value = S /> */
  bool CheckElement(const char*) const ; /** Returns TRUE if the element is present in the XML document */
      
 protected:
  TiXmlDocument *_fpXmlDoc;
};

#endif
