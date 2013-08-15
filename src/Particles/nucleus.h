/**
   @file    nucleus.h
   @author  Eric Armengaud, armengau@in2p3.fr & Nils Nierstneh�fer nils@physik.uni-wuppertal.de & Joerg Kulbartz jkulbart@mail.desy.de
   @brief   Class describing propagating protons and neutrons
*/

/*
changed by Nils Nierstenhoefer to propagate nuclei in 09/2009
*/

#ifndef _NUCLEUS_H_
#define _NUCLEUS_H_

#include "particlepropa.h"

#include "neutrino.h"
#include "photon.h"
#include "list1dphotons.h"

#include "env1d.h"
#include "units.h"
#include <fstream>
//#include <string.h>
#include <math.h>
#include "CLHEP/Random/RandFlat.h"
#include "interactiondata.h"
#include "basicpinteractions.h"
#include "sophiainteractions.h"
#include "TabulatedTALYSAveragedCrossSection.h"
//#include "universe.h"
#include "continuoussources.h"
#include "nucleusdb.h"
#include <TabulatedTALYSMeanFreePath_CMB.h>
#include <TabulatedTALYSMeanFreePath_IRB.h>


#include "photonspectrum.h"

#ifdef HAVE_TROOT_H
#include "TF1.h"
#endif

using namespace std;

/**
   @class TNucleus
   @brief Protons and neutrons as an implementation of the TParticlePropa class

   All interactions with the relevant backgrounds are implemented. Secondary particles are generated if required by the configuration file.

*/

class TNucleus : public TParticlePropa {

 public:
  TNucleus(int aMassNumber, int aChargeNumber) ;
  TNucleus(TUniverse*, TList1DPhotons*
		   ) ;
  /**< Standard constructor from a source. Generates a nucleus according to the simulation parameters. */
  TNucleus(TUniverse* aUniv, 
	   TVector3D aPosition, 
	   TVector3D aMomentum, 
	   TVector3D aInitPosition, 
	   TVector3D aInitMomentum,
	   double aTime, 
	   int aMassNumber, 
	   int aChargeNumber, 
	   TList1DPhotons* apList1DPhotons,
	   int aInitType
	   ) ;
  /**< Constructor used when a nucleus is generated as a secondary. */
  ~TNucleus() {} ;
  QUEUE<TParticle*>* Propagate( TBasicParam* aBasic) ;
  void PairProduction( double lTimeStep) ;
  /**< Pair production is handled as a continuous energy loss mechanism, using pre-computed energy loss tables and taking into account simple redshift effects. */
  void PionInteraction() ;
  /**< Stand-alone pion interaction routine, using pre-computed interaction rate tables. Checks if there is pion production during a given time step. */
  void PionProduction(int) ;
  /**< Routine called by PionInteraction. Computes the energy loss of a nucleus when pion production takes place, from energy loss tables, and converts proton <-> neutron if necessary. */
  //  void NeutronDecay() ;
  /**< Checks for neutron decay. If required, generates secondary particles associated to the decay. */
  //void SophiaInteraction() ;
  /**< Pion interaction routine for protons and neutrons on the CMB, using the SOPHIA event generator. */
  //void SophiaIRInteraction() ;
  /**< Pion interaction routine for protons and neutrons on the IRB, using the SOPHIA event generator. */
  void SophiaPionProd( const int, int aNucleonChargeNumber, double EnergyPerNucleon) ;
  /**< Interface with the event generator SOPHIA (only for protons and photons) (written in Fortran). Generates secondaries if required. */
  
  //Changed by NN
  inline double Mass();
  /**< Particle mass has to be calculated more carefuly for nuclei beause of the massdeficit which appears due to the binding energy. */
  void GetTalysPhotoDisintegrationProbability(); 
  /**< Calulates the probability for a spallation process  in a photonfield (for example in CMB)(only for nuclei with A>=12 where A is the mass number). */
  void GetTalysPhotoDisintegrationProducts();   
  /**< Generates secondaries for a spallation process  in a photonfield - if required (for example in CMB)(only for nuclei with A>=12 where A is the mass number). */
/*   double GetNuclearDecayProbability( int aMassNumber, int aChargeNumber, double dt); */
  void TALYSNucleiPhotoDisintegration(TBasicParam* aBasic);
  /**< Calculates the probability of the nucleus to decay.*/
  void GetNuclearDecayProducts();
  /**< Generates the products of a decay if required.*/
 void BetaMinusDecay();
    /**<Handels beta- Decay.*/
  void BetaPlusDecay();
  /**<Handels beta+ Decay.*/
  void AlphaDecay();
  /**<Unsurprisingly handels alpha decay.*/
  void ProtonDripping();   
  /**<Handels proton dripping decay.*/
  void NeutronDripping();
  /**<Handels neutron dripping decay. */
  //  void NeutronDecay();
  //void TRY_SophiaInteraction( int caseFlag);
  //Changed by JK
  //void GetNuclearDecayHalfLife();
  /**< Returns the half life of a nucleus in it's rest frame. Can be used to identify "stable" nuclei which will be propagated while "unstable" ones might be converted directly into the following "stable" nucleus before the next propagation step. */
  double DecayGammaTau();
  //Up to now simple copy of TNucleus::NeutronDecay with some modifications  
  void Interact(){throw TCrpErr("Interact() gerufen, obwohl das nicht passieren duerfte!"); };
  void CreatePionTable(void);
  /**< At each step, change energy and momentum norm. This takes into account all interactions, and generates secondaries if necessary. */
  
  double GetHalfLife(int aMassNumber,int aChargeNumber);
  //Returns the Halflife for a given isotope. 
  
  double PionMFP( int caseFlag);
  double PionNeutronMFP( int caseFlag);
  double PionProtonMFP( int caseFlag);
  double PairProdLossLength();

  //int _fcounter;
  
  double SibyllInteraction();
  double CalcIRLossRate(double lRefEnergy, TPhotonSpectrum infra);
  double functs_int_ir(double* eps_ln, double* param);
  void SibyllPProd();
  double RandPowerLaw(double index, double min, double max, double toto);
  
 private :
  TList1DPhotons* _fpList1DPhotons ;
  
};

#endif
