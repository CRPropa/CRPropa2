/**
	@file   units.h
	@author Tristan Beau, beau@in2p3.fr
	@brief  Usefull add ons to the CLHEP SystemOfUnit.
*/

#ifndef _UNITS_H_
#define _UNITS_H_

#include "CLHEP/Units/SystemOfUnits.h"
#include "CLHEP/Units/PhysicalConstants.h"
using namespace CLHEP;

// ----------- Time   ----------------
static const double minute             = 60 * second ;
static const double hour               = 60 * minute ; 
static const double day                = 24 * hour ;
static const double year               = 365.25 * day ;
static const double kiloyear           = 1.e3 * year; 
static const double megayear           = 1.e6 * year ;
static const double gigayear           = 1.e9 * year ;

static const double yr                 = year;
static const double kyr                = kiloyear;
static const double Myr                = megayear ;
static const double Gyr                = gigayear ;

// ----------- Length ----------------
static const double lightyear          = c_light * yr ; 
static const double inv_c_light        = 1. / c_light ;

static const double kiloparsec         = 1.e3 * pc; 
static const double megaparsec         = 1.e6 * pc;
static const double gigaparsec         = 1.e9 * pc;

static const double kpc                = kiloparsec ;
static const double Mpc                = megaparsec ;
static const double Gpc                = gigaparsec ;

static const double inverse_megaparsec = 1. / megaparsec;
static const double inverse_kiloparsec = 1. / kiloparsec;
static const double inverse_gigaparsec = 1. / gigaparsec;

static const double invkpc             = inverse_kiloparsec;
static const double invMpc             = inverse_megaparsec;
static const double invGpc             = inverse_gigaparsec;

// ----------- Energy ----------------
static const double exaelectronvolt    = 1.e9 * GeV   ;
static const double zetaelectronvolt   = 1.e12 * GeV  ;

static const double EeV                = exaelectronvolt ;
static const double ZeV                = zetaelectronvolt ; 

static const double invEeV             = 1. / EeV;
static const double invGeV             = 1. / GeV;

// -------- Magnetic fields ---------
static const double microgauss         = 1.e-6 * gauss ;
static const double nanogauss          = 1.e-9 * gauss ;

static const double muG                = microgauss ;
static const double nG                 = nanogauss ;

#endif
