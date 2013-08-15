/**
   @file   continuoussources.h
   @author Eric Armengaud, armengau@in2p3.fr
   @brief  Class describing continuous sources
*/

#ifndef _CONTINUOUSSOURCES_H
#define _CONTINUOUSSOURCES_H

#include "sources.h"
#include "xmlparam.h"
#include "crp_err.h"

#include "fitsio.h"

#include "fits_err.h"

/**
   @class TContinuousSources
   @brief Class describing continuous sources: a 1-D or 3-D source density must be specified. During the simulation, each trajectory will start from a point choosen at random from this density.
*/

class TContinuousSources : public TSources, public TXmlParam {
 public:
  TContinuousSources(const char*) ;
	~TContinuousSources() ;

  TVector3D GetPosition() const { return _fpSourceDensity->getSourcePosition(); }
  /**< Returns a position drawn at random from the source density. */

};

#endif
