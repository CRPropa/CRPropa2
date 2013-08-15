
/*
 *  typedclass.h
 *  CRPropa_head
 *
 *  Created by Tristan Beau (beau@in2p3.fr) on 17/11/2005.
 *  Copyright 2005 APC-IN2P3. All rights reserved.
 *
 *  Nils Nierstenhoefer (NN) added PARTICLE_NUCLEUS in 09/2009
 *
 */

#ifndef _TYPEDCLASS_H_
#define _TYPEDCLASS_H_


enum UniverseType {
  UNIVERSE_ENV1D = 0x101,
  UNIVERSE_GALACTIC,
  UNIVERSE_LARGESCALE
} ;

enum SourceType {
  SOURCE_CONTINUOUS = 0x201,
  SOURCE_DISCRETE
} ;

enum OutputType {
  OUTPUT_ASCII = 0x301,
  OUTPUT_FITS,
  OUTPUT_ROOT,
  OUTPUT_NO
} ;

enum ParticleType {
  PARTICLE_NEUTRINO = 0x401,
  //TODO delete : Should be obselete as a nucleon is a nucleus, too.
  PARTICLE_NUCLEON,
  PARTICLE_NUCLEUS,//Changed by NN
  PARTICLE_PHOTON
};

enum ObserverType {
  OBSERVER_LARGESPHERE = 0x501,
  OBSERVER_POINT,
  OBSERVER_SMALLSPHERE,
  OBSERVER_NO
};

enum MagFieldType {
  MAGFIELD_1D = 0x601,
  MAGFIELD_KOLMOGOROFF,
  MAGFIELD_LSS,
  MAGFIELD_NO,
  MAGFIELD_UNIFORM,
  MAGFIELD_GRID
} ;

enum InteractionType {
  INTERACTION_BASIC = 0x701,
  INTERACTION_SOPHIA,
  INTERACTION_PHOTON,
  INTERACTION_NUCLEUS,// Changed by NN
  INTERACTION_NO
} ;

enum InfraredType {
  IR_UNIFORM_PRIMACK = 0x801,
  IR_UNIFORM_HIGH,
  IR_UNIFORM_LOW,
  SHELL,
  THREED
};

enum PhotonBG_Type {
  PhotonBG_CMB = 0x901,
  PhotonBG_IR,
  PhotonBG_Opt
};

enum GasType {
  CLUSTER = 0x901,
  GRIDGAS,
  LSSGAS
};

enum PD_TabDataType{undefined,
		    CMB_MFPTab,
		    IRB_MFPTab,
		    AveragedXSTab,
		    XSTab};

/**
   @class TTypedClass
   @brief Basic class to give a type to objects that inherit from virtual classes

   This allows to use the functions Type() and SetType() for all classes inheriting from a TTypedClass. 
   Type() returns an int encoding the type of object, following codes that are given in the enumerations of
   the form ****Type (eg. SourceType)
*/

class TTypedClass {
public:
	TTypedClass(int aType=0) { SetType(aType); }
	~TTypedClass() { }

	int Type() const { return GetType(); }
	int GetType() const { return _fType;}
protected:
	int SetType(int aType) { return (_fType=aType) ; }
private:
	int _fType;
} ;

#endif
