
/**
   @file   largesphereobservers.cc
   @author Eric Armengaud, armengau@in2p3.fr
   @brief  Implementation of the TLargeSphereObservers class. See the .h file
*/

#include "largesphereobservers.h"

TLargeSphereObservers::TLargeSphereObservers(const char* aFileName) : TXmlParam(aFileName) {
  SetType(OBSERVER_LARGESPHERE);
  TiXmlElement* lpXmlObs = XmlExtract().GetElement("Observers") ;
  if (!lpXmlObs) throw TXmlErr("There should be an Observer in xml file") ;
  string lDum = lpXmlObs->Attribute("type") ;
  if (lDum != "Spheres around Source") throw TXmlErr("Bad observer type.") ;
  std::cout<< "Obs type : " << this->Type() << std::endl;
  // 1 ) Number of observers
  TiXmlElement* lpXmlNumber = lpXmlObs->FirstChildElement("Number") ;
  int lNb ;
  lDum = lpXmlNumber->Attribute("value",&lNb) ;
  if ( lNb <= 0 ) {
    throw TCrpErr("Error getting number of observers");
  }
  _fNb = lNb ;

  // 2 ) Positions and radii
  TiXmlElement *lpXmlCoordX, *lpXmlCoordY, *lpXmlCoordZ, *lpXmlRad ;
  double lX,lY,lZ,lRadius ;
  TVector3D lPos ;
  TiXmlElement* lpSphereObs = lpXmlObs->FirstChildElement("Sphere") ;
  if (!lpSphereObs) {
    throw TXmlErr("Error : no sphere specified in xml file");
  }
  for (int j=0; j<_fNb; j++) {
    lpXmlCoordX = lpSphereObs->FirstChildElement("CoordX_Mpc") ;
    if (!lpXmlCoordX) throw TXmlErr("Error : no CoordX for sphere in xml file"); 
    lDum = lpXmlCoordX->Attribute("value",&lX) ;
    lpXmlCoordY = lpSphereObs->FirstChildElement("CoordY_Mpc") ;
    if (!lpXmlCoordY) throw TXmlErr("Error : no CoordY for sphere in xml file");
    lDum = lpXmlCoordY->Attribute("value",&lY) ;
    lpXmlCoordZ = lpSphereObs->FirstChildElement("CoordZ_Mpc") ;
    if (!lpXmlCoordZ) throw TXmlErr("Error : no CoordZ for sphere in xml file");
    lDum = lpXmlCoordZ->Attribute("value",&lZ) ;
    lPos.set(lX*Mpc,lY*Mpc,lZ*Mpc) ;
    _fPositions.push_back(lPos) ;

    lpXmlRad = lpSphereObs->FirstChildElement("Radius_Mpc") ;
    if (!lpXmlRad) throw TXmlErr("Error : no Radius for sphere."); 
    lDum = lpXmlRad->Attribute("value",&lRadius) ;
    if ( lRadius <= 0 ) throw TXmlErr("Error getting radius of observers");
    _fRadii.push_back(lRadius*Mpc) ;

    lpSphereObs = lpSphereObs->NextSiblingElement("Sphere") ;
    if ( (j!=_fNb-1 && ! lpSphereObs) || (j==_fNb-1 && lpSphereObs) )
      throw TXmlErr("Error parsing sphere observers with the specified number"); 
  }

}
