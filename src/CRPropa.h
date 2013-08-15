/**
   @file   CRPropa.h
   @author Tristan Beau, beau@in2p3.fr
   @brief  Header file of the main. Define includes and error enum. 
*/

#ifndef _CRPROPA_H_
#define _CRPROPA_H_
//#define DEBUG_OUTPUT

//C includes
#include <stdlib.h>

//C++ includes
#include <iostream>
#include <fstream>
#include <limits>

// CLHEP includes
//#include "CLHEP/HepMC/GenParticle.h"
//#include "CLHEP/HepPDT/DefaultConfig.hh"
//#include "CLHEP/HepPDT/TableBuilder.hh"
//#include "CLHEP/HepPDT/ParticleDataTableT.hh"
#include "CLHEP/Random/RandGeneral.h"
#include "CLHEP/Random/RandFlat.h"

// CFITSIO includes
#include "fitsio.h"

//Local files to include
#include "defqueue.h"
#include "units.h"
#include "everything.h"
#include "outputdata.h"
#include "particle.h"
#include "env1d.h"
#include "largescalestructure.h"
#include "galacticstructure.h"
#include "sources.h"
#include "discretesources.h"
#include "observers.h"
#include "pointobserver.h"
#include "vector3d.h"
#include "nucleus.h"
#include "list1dphotons.h"
#include "interactiondata.h"
#include "basicpinteractions.h"
#include "crp_err.h"
#include "xml_err.h"
#include "sysdep.h"

using namespace std;
//using namespace HepMC;
//using namespace HepPDT;

#ifdef NILS
//For Nils' test plots (begin)
TCanvas* MassChargePathCan = new TCanvas("MassChargePathCan","MassChargePathCan",1);
TH2F* MassChargePathTH2F   = new TH2F("MassChargePathCanTH2F","MassChargePathCanTH2F",
				      27 ,0., 27.,
				      57, 0., 57. );
//TCanvas* MassChargePathAddOnCan = new TCanvas("MassChargePathAddOnCan","MassChargePathAddOnCan",1);
/*TH2F* MassChargePathAddOnTH2F   = new TH2F("MassChargePathAddOnCanTH2F","MassChargePathAddOnCanTH2F",
				      27*500 ,0., 27.,
				      57*500, 0., 57. );*/

TCanvas* MassChargeVSChannelCan = new TCanvas("MassChargeVSChannelCan","MassChargeVSChannelCan",1);
TH2F* MassChargeVSChannelTH2F = new TH2F("MassChargeVSChannelCanTH2F","MassChargeVSChannelCanTH2F",
					 15*5*2 ,5, 20,
					 9, 0., 9. );
TCanvas* MFPvsEnCan = new TCanvas("MFPvsEnCan","MFPvsEnCan",1);
TH2F* MFPvsEnTH2F = new TH2F("MFPvsEnCanTH2F","MFPvsEnCanTH2F",
			     35 , 6, 14,
			     200, -5., 9. );

TCanvas* MFPvsEn_ExclCan = new TCanvas("MFPvsEn_ExclCan","MFPvsEn_ExclCan",1);
TH2F* MFPvsEn_ExclTH2F = new TH2F("MFPvsEn_ExclCanTH2F","MFPvsEn_ExclCanTH2F",
			     35 , 6, 14,
			     200, -5., 9. );

int NilsEvCounter=0.;
//For Nils' test plots (end)
#endif


#endif
