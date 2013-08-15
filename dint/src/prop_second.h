
#ifndef _PROPSECOND_H_
#define _PROPSECOND_H_

#include "cvector.h"

#include <stdio.h>
#include <math.h>
#include "rate.h"
#include "const.h"
#include "spectrum.h"
#include "cvector.h"
#include "load.h"
#include "prepare.h"
#include "sync.h"
#include "inject.h"
#include "background.h"
#include "fold.h"
#include "advance.h"
#include "final.h"
#include "utilities.h"
#include <string.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

/*--------------------------------------------------------------------------
  prop_second
  Main routine for the dint package
  E. Armengaud - Jan 2006
--------------------------------------------------------------------------*/

/*
  Arguments :

  dist_observer : in Mpc
  pB_field : in Gauss
  pEnergy, pEnergyWidth : in ELECTRON_MASS*eV
  apInjectionSpectrum : INPUT mean number of particles per energy bin
  pSpectrum : OUTPUT mean number of particles per energy bin
  aDirTables : directory for cascade tables
  Photon background flags :
       aIRFlag = 0 (High IR), 1 (Low IR) or 2 (Primack IR)
       aZmax_IR = max redshift of the infrared background
       aRadioFlag = 0 (High Radio), 1 (Low Radio), 2 (Obs Radio) or 3 (Null Radio)
  Cosmological parameters : H0, Omega_M/Lambda. The unit of H0 is km/s/Mpc
  aCutcascade_Magfield : flag to "cut" of the e+/- cascade by the magnetic deflections.
*/

/*--------------------------------------------------------------------------
 The routine prop_second, used in CRPropa, propagates an electromagnetic cascade
 in the extragalactic medium over a given distance (redshifts taken into account)

 Propagation takes place in redshift space. The energy unit used in dint is m_e (~511 keV).
 The distance unit is in cm.
 The magnetic field perpendicular to the trajectory must be specified. It can be 
 inhomogeneous, allowing to take into account probable B field concentrations inside the
 clusters. 3 models of homogeneous cosmic IR background and 3 models of cosmic radio
 background are implemented.

 The redshift evolution of radio background is consistent with its model (see background.cpp).
 The redshift evolution of IR background is the same as for CMB until aZmax_IR.

 The input spectrum must have the same format as the output : 
 The output is a spectrum, computed for 170 values of energy ranging from 10^7 to 10^24 eV.

 If the aCutcascade_Magfield flag is set : at each step, we define a critical energy E_c 
 based on the comparison between Larmor, synchrotron and ICS time scales.
 Only the electrons with  r_Larmor > aCutcascade_Magfield) * min(t_synchrotron,t_ICS) are kept:
 this mimics the loss of the cascade due to the deflections (useful for the study of a "point" 
 source; allows to test roughly the effect of B fields on the 1D approximation for the cascade).
--------------------------------------------------------------------------*/

void prop_second(const double dist_observer, 
		 //const double InjEnergy,
		 //const PARTICLE part, 
		 //const double HInjEnergy, const double deltaE_hadron,
		 const dCVector* pB_field,
		 const dCVector* pEnergy,
		 const dCVector* pEnergyWidth,
		 Spectrum* apInjectionSpectrum,
		 Spectrum* pSpectrum, string aDirTables,
		 const int aIRFlag, const double aZmax_IR, const int aRadioFlag,
		 const double aH0, const double aOmegaM, const double aOmegaLambda,
		 const double aCutcascade_Magfield);

void BuildRedshiftTable(double aH0, double aOmegaM, double aOmegaLambda, 
			dCVector* pRedshiftArray, dCVector* pDistanceArray);

double getRedshift(dCVector RedshiftArray, dCVector DistanceArray, double distance) ;
double getDistance(dCVector RedshiftArray, dCVector DistanceArray, double redshift) ;

#endif
