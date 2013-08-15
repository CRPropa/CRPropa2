
/**
   @file   smallsphereobservers.cc
   @author Eric Armengaud, armengau@in2p3.fr
   @brief  Implementation of the TSmallSphereObservers class. See the .h file
*/

#include "smallsphereobservers.h"

TSmallSphereObservers::TSmallSphereObservers(const char* aFileName) : TXmlParam(aFileName) {
SetType(OBSERVER_SMALLSPHERE);
  TiXmlElement* lpXmlObs = XmlExtract().GetElement("Observers") ;
  if (!lpXmlObs) throw TCrpErr("There should be an Observer in xml file") ;
  string lDum ;

  // 1 ) Number of observers
  TiXmlElement* lpXmlNumber = lpXmlObs->FirstChildElement("Number") ;
  int lNb ;
  lDum = lpXmlNumber->Attribute("value",&lNb) ;
  if ( lNb <= 0 ) {
    throw TCrpErr("Error getting number of observers");
  }
  _fNb = lNb ;

  // 2 ) Radius
  TiXmlElement* lpXmlRad = lpXmlObs->FirstChildElement("Radius_Mpc") ;
  double lRadius ;
  lDum = lpXmlRad->Attribute("value",&lRadius) ;
  if ( lRadius <= 0 ) {
    throw TCrpErr("Error getting radius of observers"); 
  }
  _fRadius = lRadius * Mpc ;

  // 3 ) Positions
  TiXmlElement *lpXmlCoordX, *lpXmlCoordY, *lpXmlCoordZ ;
  double lX,lY,lZ ;
  TVector3D lPosition ;
  TiXmlElement* lpSphereObs = lpXmlObs->FirstChildElement("SphereObserver") ;
  if (!lpSphereObs) throw TCrpErr("Error : no sphere observer specified in xml file");
  for (int j=0; j<_fNb; j++) {
    lpXmlCoordX = lpSphereObs->FirstChildElement("CoordX_Mpc") ;
    if (!lpXmlCoordX) throw TCrpErr("Error : no CoordX for sphere observer in xml file"); 
    lDum = lpXmlCoordX->Attribute("value",&lX) ;
    lpXmlCoordY = lpSphereObs->FirstChildElement("CoordY_Mpc") ;
    if (!lpXmlCoordY) throw TCrpErr("Error : no CoordY for sphere observer in xml file");
    lDum = lpXmlCoordY->Attribute("value",&lY) ;
    lpXmlCoordZ = lpSphereObs->FirstChildElement("CoordZ_Mpc") ;
    if (!lpXmlCoordZ) throw TCrpErr("Error : no CoordZ for sphere observer in xml file"); 
    lDum = lpXmlCoordZ->Attribute("value",&lZ) ;
    lPosition.set(lX*Mpc,lY*Mpc,lZ*Mpc) ;
    _fPositions.push_back(lPosition) ;
#ifdef DEBUG_OUTPUT 
    std::cout<< "Read in small sphere observer # "<<j<<" at ("
	     <<lPosition.x()/Mpc<<","
	     <<lPosition.y()/Mpc<<","
	     <<lPosition.z()/Mpc<<")"
	     <<" / F.Y.I. unit conversion via Mpc="<<Mpc
	     <<std::endl;
#endif
    lpSphereObs = lpSphereObs->NextSiblingElement("SphereObserver") ;
    if ( (j!=_fNb-1 && ! lpSphereObs) || (j==_fNb-1 && lpSphereObs) ) {
      throw TCrpErr("Error parsing point sources with the specified number in xml file"); 
    }
  }

}
