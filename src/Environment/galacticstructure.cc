
/**
   @file    galacticstructure.cc
   @author  Eric Armengaud, armengau@in2p3.fr
   @brief   Implementation of the TGalacticStructure class. See the .h file
*/

#include "galacticstructure.h"

TGalacticStructure::TGalacticStructure(const char* aFileName) : TXmlParam(aFileName) {
	SetType(UNIVERSE_GALACTIC);
	throw TCrpErr("Galactic environment not implemented.");
}

