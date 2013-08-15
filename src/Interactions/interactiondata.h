/**
   @file    interactiondata.h
   @author  Eric Armengaud, armengau@in2p3.fr
   @brief   Virtual class describing the parameters for any interaction
*/

#ifndef _INTERACTIONDATA_H_
#define _INTERACTIONDATA_H_

#include <iostream>
#include <vector>
#include "units.h"
#include "CLHEP/Random/RandGeneral.h"
#include "TabulatedTALYSY.h"
#include "typedclass.h"
//#include "NucleusDB.h"
#include "prop_second.h" // pour NUM_MAIN_BINS

#define NBINS_PAIRPROD 70 // Used in sofia and photon interactions, for pair prod
#define EPMIN_PAIRPROD 15
#define EPMAX_PAIRPROD 22

//#ifdef HAVE_TROOT_H
#include <TGraph.h>
//#endif
#include "vector3d.h"

using namespace std;

/**
   @class TInteractionData
   @brief Virtual class for interaction parameters

*/

class TInteractionData : public TTypedClass {
 public:
  
  virtual ~TInteractionData() {}
  
  double InteractionTimeStep() const { return _fInteractionTimeStep; }
  /**< This is the maximum time step allowed by interactions. It is taken into account during propagation. */

  virtual string ShowerTableDir() const { return 0; }

  virtual double L_pair(int) const { return 0; }
  virtual double E_part(int) const { return 0; }
  virtual double dEtabbin() const { return 0; }
  virtual double inc_p(int) const { return 0; }
  virtual double inc_n(int) const { return 0; }
  virtual double Eminp(int) const { return 0; }
  virtual double Eminn(int) const { return 0; }
  virtual double lossp(int, int) const { return 0; }
  virtual double lossn(int, int) const { return 0; }
  virtual double Lpp(int, int) const { return 0; }
  virtual double Lnn(int, int) const { return 0; }
  virtual double lossp_tot(int) const { return 0; }
  virtual double lossn_tot(int) const { return 0; }

  virtual double E_pionprod(int) const { return 0; }
  virtual double E_IRpionprod(int) const { return 0; }
  virtual double LossRateProton(int) const { return 0; }
  virtual double LossRateNeutron(int) const { return 0; }
  virtual double IRLossRateProton(int) const { return 0; }
  virtual double IRLossRateNeutron(int) const { return 0; }
  virtual double GetPPEps() const {return .01;}

  virtual double dEtabPion() const { return 0; }
  virtual double dEtabIRPion() const { return 0; }

  virtual bool SecPairProdPhotonFlag() const { return 0; }
  virtual double SecPairProdPhotonProba() const { return 0; }
  virtual bool PairProdFlag() const { return 0; }
  virtual bool PionProdFlag() const { return 0; }
  virtual bool DecayFlag() const {return 0;}
  virtual bool KelnerPairProdFlag() const { return 0; }
  virtual bool RedshiftFlag() const { return 0; }
  virtual bool RedshiftEnergyLossFlag() const { return 0; }
  virtual bool IRPionProdFlag() const { return 0; }
  virtual bool PhotodisintegrationFlag() const { return 0; }
  virtual bool SecPhotonFlag() const { return 0; }
  virtual bool SecNuFlag() const { return 0; }
  virtual double CutcascadeFlag() const { return 0; }
  virtual double PairProdSpec(double, int) const { return 0; }

  virtual RandGeneral* RandDistriNeutronDecay() const { return 0; }

  virtual vector <unsigned long int> GetNucleusTransition(unsigned long int, double, double , double, double){}
  /* 
  virtual unsigned long int GetExclusiveChannel(unsigned long int,
						 unsigned long int, 
						 double,
						 double,
						 int,
						 int,
						 double,
						 double){}*/

  //#ifdef HAVE_TROOT_H
   //TODO : DELETE This Memeber : RootPlotMeanFreePathVsEnergy
  virtual void RootPlotMeanFreePathVsEnergy(unsigned long int,
					    TVector3D,
					    double,
					    double,
					    double,
					    double,
					    int,
					    double){}

  /* virtual void RootFileWithAllMFPs(TabulatedTALYSY* ,
				   double,
				   double,
				   double,
				   double ){}
  virtual TGraph* RootPlotYVsX(TabulatedTALYSY* ,
			       unsigned long int,
			       long unsigned int,
			       double,
			       double,
			       double,
				double,
			       int){}
  virtual TGraph* RootPlotExclSumYVsX(TabulatedTALYSY* ,
				      unsigned long int,
				      double,
				      double,
				      double,
				      double){}*/
  //#endif

   virtual unsigned long int GetExclusivePDChannelDirectly(unsigned long int InitNucleusId,
							   TVector3D Position,
							   double Energy,
							   double NucleusMass,
							   double z,
							   double TimeStep){}
   
   virtual double GetPDTimeStep(unsigned long int,
				TVector3D,
				double,
				double,
				double){}
   virtual std::vector<TabulatedTALYSY*> GetPointerToPhotoDisTables(){}
   
   virtual double XSecPPProton(int) const { return 0; }
   virtual double XSecPPProton(double)  { return 0; }
   virtual double IRXsecInt_Proton(int,int) const { return 0; }
   virtual double IRXsecInt_Neutron(int,int) const { return 0; }
   virtual double E_IR(int) const { return 0; }
   virtual double S_IR(int) const { return 0; }
   virtual double dEtabXsecIR() const { return 0; }
   virtual double dStabXsecIR() const { return 0; }
   virtual bool PProdFlag()    const { return 0; }
   

 protected:
  double _fInteractionTimeStep;
  //const TNucleusDB* _fpNucleusDB;
}; 

#endif
