/**
   @file   vector3d.h
   @author Tristan Beau, beau@in2p3.fr
   @brief  Class of 3d vectors, currently inheritated from CLHEP vector class
*/

#ifndef _VECTOR3D_H_
#define _VECTOR3D_H_

#include "CLHEP/Vector/ThreeVector.h"
using namespace CLHEP;

typedef Hep3Vector TVector3D ;

/*inline std::ostream & operator<<(std::ostream& aS, const TVector3D& aV ) {
  aS << aV.x() << "\t" << aV.y() << "\t" << aV.z() ;
  return aS;
}
*/

#endif
