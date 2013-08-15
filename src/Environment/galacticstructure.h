/**
   @file    galacticstructure.h
   @author  Eric Armengaud, armengau@in2p3.fr
   @brief   Class describing a galactic-scale environment
*/

#ifndef _GALACTICSTRUCTURE_H_
#define _GALACTICSTRUCTURE_H_

#include "universe.h"
#include "xmlparam.h"

/**
   @class TGalacticStructure
   @brief Galactic environment, not yet implemented

*/

class TGalacticStructure : public TUniverse, public TXmlParam {
 public:
	TGalacticStructure(const char*) ;
  ~TGalacticStructure() {} ;
 // int Type() const { return UNIVERSE_GALACTIC; } ;

};

#endif
