#/**
   @class  IRBzEvolutionFactory

   @author  Nils Nierstenhoefer nierstenhoefer@physik.uni-wuppertal.de

   @brief  Singleton+factory class to provide the scaling factor which rescales the mean free path \f$\lambda \f$ of a particle interaction in the IRB at redshift z=0 to a mean free path at redhsift z'>0. 

   It is assumed that only one IRB redshift model is needed for one instance of CRPropa. Thus, IRBzEvolutionFactory is organized as a singelton. 
   In this way it can easily be called from all parts of CRPropa. But, the first call has to be performed with GetInstance(std::string XMLFileName) 
   which is currently done in everything.cc. After creation of the instance of the singelton, one can alternatively use GetInstance()->GetScalingFactor(redshift) 
   to get the scaling factor for a given redshift z'.  
   
   The main task of  this class is to create an instance of a concrete evolution model derived from TIRBzEvolutionModel (such a pattern/behaviour would presumably 
   be called "a factory" amongst programmers). Currently, there are three models to choose from:
   - Kneiske2004_BestFit
   - CMBLikeScaling
   - NoIRBScaling

   Note, that TIRBzEvolutionFactory::GetScalingFactor(double redshift) simply calls the TIRBzEvolutionModel::GetScalingFactor(double redshift) of the concrete model which was choosen.

   TIRBzEvolutionFactory and TIRBzEvolutionModel inherit from TXmlParam. Thus, their instances configure themselves automatically according to the specifications in the xml file. 

   The default model is Kneiske2004_BestFit (this is consistent as currently this Kneiske IRB parametrization is is the default IRB). 
*/

#ifndef _IRBzEvolutionFactory_H_
#define _IRBzEvolutionFactory_H_
#include "xmlparam.h"
#include "xmlextract.h"
#include <TGraph.h>
#include <string>
#include <TIRBzEvolutionModel.h>
#include <NoIRBScaling.h> 
#include <Kneiske2004_BestFit.h> 
#include <CMBLikeScaling.h> 
#include <UserDefinedScaling.h>
#include <TLegend.h>
#include <TCanvas.h>

class IRBzEvolutionFactory : public TXmlParam{
 public:
  static IRBzEvolutionFactory* GetInstance(std::string XMLFileName);
  static IRBzEvolutionFactory* GetInstance();
  double GetScalingFactor(double redshift){
    return _fModelInstance->GetScalingFactor(redshift);
  };
  void PlotAllModels();
  
 private:
  static IRBzEvolutionFactory* _fFactoryInstance;
  static TIRBzEvolutionModel* _fModelInstance;
  static std::string _fFileName;
  IRBzEvolutionFactory(std::string XMLFileName): TXmlParam(XMLFileName.c_str()){};            
  ~IRBzEvolutionFactory();                
  IRBzEvolutionFactory( const IRBzEvolutionFactory& );              
  IRBzEvolutionFactory* operator = (const IRBzEvolutionFactory &); 
  static TIRBzEvolutionModel* GetConcreteIRBzEvolutionModel(std::string XMLFileName);
};
#endif
