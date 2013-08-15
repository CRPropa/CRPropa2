/**
   @class TIRBzEvolutionModel

   @author  Nils Nierstenhoefer nierstenhoefer@physik.uni-wuppertal.de

   @brief Abstract class to model IRB redshift evolution by a simple rescaling of the mean free paths of particle interactions. 

   As outlined in the CPropa paper for version 1.4 (page 12 equation 13) one can scale the interaction length \f$\lambda \f$ of a particle interaction in the CMB at redshift z=0 to yield \f$\lambda' \f$ at redshift z'
 
   \f$\lambda'[E]^{-1} = S(z) \times \lambda[(1+z')E]^{-1}\f$

with

   \f$ S(z)=(1+z)^3\f$.

   Classes can be derived from TIRBzEvolutionModel that provide a similar scaling function S(z) for interactions with IRB photons. In this case the scaling factor is defined as
   
\f$ S(z) = (1+z)^3 \frac {\int_{0}^{\infty} n_{IRB}(\epsilon, z) d \epsilon} {\int_{0}^{\infty} n_{IRB}(\epsilon, z=0) d \epsilon} \f$

where \f$ n(\epsilon, z) \f$, \f$ \epsilon \f$ are the photon number denistiy [GeV^2], photon energy repectively - both in comoving units. 

   Note, instances of classes derived from TIRBzEvolutionModel (e.g. Kneiske2004_BestFit) should not be called directly. Use IRBzEvolutionFactory if you want to access a concrete IRB scaling model. 
*/

#ifndef _TIRBzEvolutionModel_H_
#define _TIRBzEvolutionModel_H_
#include "xmlparam.h"
#include "xmlextract.h"
#include <TGraph.h>
#include <PhotonBackground.h>

class TIRBzEvolutionModel: public TXmlParam{
 public:
  virtual double GetScalingFactor(double redshift)=0;
  //The only way to create a IRB evolution model is via factory which is a friend class!
  TIRBzEvolutionModel(std::string XMLFileName);   
  //  TIRBzEvolutionModel( const TIRBzEvolutionModel& );             
  //TIRBzEvolutionModel & operator = (const TIRBzEvolutionModel &); 
  TGraph PlotModelTGraph(double zMin, double zMax, int NPoints);
  virtual ~TIRBzEvolutionModel() {}

 protected:
  double _fMaxRedshift;
};

#endif
