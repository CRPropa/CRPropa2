/**
   @file    everything.h
   @author  Tristan Beau, beau@in2p3.fr
   @brief   Class which will contain every usefull info for simu
*/

#ifndef __EVERYTHING_H__
#define __EVERYTHING_H__

#include "universe.h"
#include "env1d.h"
#include "galacticstructure.h"
#include "largescalestructure.h"
#include "basicparam.h"
//#include "IRBzEvolutionFactory.h"
/**
   @class TEveryThing
   @brief This class contains all the parameters of the simulation, namely a pointer to a TBasicParam object, and a pointer to a TUniverse.
*/

class TEveryThing {
	public:
	TEveryThing(const char*);
	~TEveryThing();
	
	TBasicParam *Basic;
	TUniverse   *Univ;
};

#endif
