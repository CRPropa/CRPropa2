
/**
   @file   pointobserver.cc
   @author Eric Armengaud, armengau@in2p3.fr
   @brief  Implementation of the TPointObserver class. See the .h file
*/

#include "pointobserver.h"

TPointObserver::TPointObserver(const char* aFileName) : TXmlParam(aFileName) {
  SetType(OBSERVER_POINT);
  TiXmlElement* lpEnv = XmlExtract().GetElement("Environment") ;
  string lEnv = lpEnv->Attribute("type") ;
  if (lEnv == "One Dimension") {
    _fNb = 1 ;
  } else throw TCrpErr("Class TPointObserver only in 1D");
  
}
