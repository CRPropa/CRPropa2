
/**
   @file   xmlextract.cc
   @author Tristan Beau, beau@in2p3.fr
   @brief  Implementation of the TXmlExtract class. See .h file for other comments.
*/

#include "xmlextract.h"

#include <fstream>
using namespace std;

string GetAttribute(const TiXmlElement* aElem, const char* aAttribute) {
  const char *lpC = aElem->Attribute(aAttribute);
  if ( ! lpC ) 
    throw TXmlErr( string("Undefined attribute : ")+string(aAttribute), false );
  
  return string ( lpC );
}

TXmlExtract::TXmlExtract(const char *aFileName) {

  _fpXmlDoc = new TiXmlDocument ( aFileName );

  if ( ! _fpXmlDoc->LoadFile() ) 
    throw TXmlErr( string("Error while loading file : ") + string(_fpXmlDoc->ErrorDesc()) );
}

int TXmlExtract::GetInt(const char* aDescription) const {
  int lReturn ;
  const char *lTest;

  TiXmlElement *lElem = _fpXmlDoc->FirstChildElement(aDescription) ;
  if ( ! lElem ) throw TXmlErr(string("Element not found: ")+string(aDescription)) ;
  lTest = lElem->Attribute("value", &lReturn);
  if ( ! lTest )  throw TXmlErr(string("Missing 'int' attribute for element: ")+string(aDescription)) ;

  return lReturn;
}

unsigned long TXmlExtract::GetULong(const char* aDescription) const {
  unsigned long lReturn ;
  const char *lTest;

  TiXmlElement *lElem = _fpXmlDoc->FirstChildElement(aDescription) ;
  if ( ! lElem ) throw TXmlErr(string("Element not found: ")+string(aDescription)) ;
  lTest = lElem->Attribute("value");
  if (lTest) lReturn = strtoul(lTest, NULL, 10) ;
  else throw TXmlErr(string("Missing 'long' attribute for element: ")+string(aDescription)) ;

  if (lReturn == ULONG_MAX) throw TXmlErr("Too long number in GetULong() !") ; 
  return lReturn;
}

double TXmlExtract::GetDouble(const char* aDescription) const {
  double lReturn ;
  const char *lTest;

  TiXmlElement *lElem = _fpXmlDoc->FirstChildElement(aDescription) ;
  if ( ! lElem ) throw TXmlErr(string("Element not found: ")+string(aDescription)) ;
  lTest = lElem->Attribute("value", &lReturn);
  if ( ! lTest ) throw TXmlErr(string("Missing 'double' attribute for element: ")+string(aDescription)) ;

  return lReturn;
}

string TXmlExtract::GetString(const char* aDescription) const {
  string lReturn ;
  const char* lVal ;

  TiXmlElement *lElem = _fpXmlDoc->FirstChildElement(aDescription) ;
  if ( ! lElem ) throw TXmlErr(string("Element not found: ")+string(aDescription)) ;
  lVal = lElem->Attribute("value");
  if ( ! lVal ) throw TXmlErr(string("Missing 'string' attribute for element: ")+string(aDescription)) ;
  lReturn = lVal ;

  return lReturn;
}

TiXmlElement* TXmlExtract::GetElement(const char* aDescription) const {
  TiXmlElement *lElem = _fpXmlDoc->FirstChildElement(aDescription) ;
  if (! lElem) throw TXmlErr( string(" Element not found: ") + string(aDescription) );

  return lElem ;
}


bool TXmlExtract::CheckElement(const char* aDescription) const {
  bool lCheck = 1 ;
  TiXmlElement *lElem = _fpXmlDoc->FirstChildElement(aDescription) ;
  if ( ! lElem ) lCheck = 0 ;

  return lCheck ;
}

TXmlExtract::~TXmlExtract() {

  delete _fpXmlDoc ;

}
