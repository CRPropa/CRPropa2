/**
   @class  UseDefinedScaling
   
   @author  Peter Schiffer peter.schiffer@desy.de
   
   @brief  In this IRB evolution model the mean free path (mfp) for particle interactions with IRB photons is scaled with a user defined factor. A scaling with (1+z,3) corresponds to a CMB-like scaling. The scaling is defined in a text file specified in the xml file


*/

#ifndef _UserDefinedScaling_H_
#define _UserDefinedScaling_H_
#include <TIRBzEvolutionModel.h>

//2. The case of CRPropa 1.4 cmb-like scaling pow(1+z,3) with MaxRedshift at which the IRB will be switched off.
class UserDefinedScaling : public TIRBzEvolutionModel{
 public:
  UserDefinedScaling(std::string XMLFileName);
  double GetScalingFactor(double redshift);
  virtual ~UserDefinedScaling() { }

 protected:
  vector<double> _fRedshiftArray ; /**< Distance = angular distance = comoving distance/(1+z) */
  vector<double> _fScalingArray ;


};
#endif
