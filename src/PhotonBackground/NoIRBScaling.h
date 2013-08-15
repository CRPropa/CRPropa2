/**
   @class  NoIRBScaling
   
   @author  Nils Nierstenhoefer nierstenhoefer@physik.uni-wuppertal.de
   
   @brief  Returns 1 for all redshifts <=> no scaling.
*/

#ifndef _NoIRBScaling_H_
#define _NoIRBScaling_H_
#include <TIRBzEvolutionModel.h>

//1. NoIRBScaling object returns =1 for all z
class NoIRBScaling : public TIRBzEvolutionModel{
 public:
  NoIRBScaling(std::string XMLFileName): TIRBzEvolutionModel::TIRBzEvolutionModel(XMLFileName) {};
  double GetScalingFactor(double redshift);

  virtual ~NoIRBScaling() { }
};

#endif
