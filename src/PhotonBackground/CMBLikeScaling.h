/**
   @class  CMBLikeScaling
   
   @author  Nils Nierstenhoefer nierstenhoefer@physik.uni-wuppertal.de
   
   @brief  In this IRB evolution model the mean free path (mfp) for particle interactions with IRB photons is scaled  with (1+z,3) as in the case of the cmb where z is the redshift. Furthermore, the scaling factor is set to zero for z>MaxRedshift where MaxRedshift can be defined via the xml settings in the IRBBackground scope. This corresponds to an instantaneous creation of the IRB at redshift MaxRedshift.

Note, that this approach has been used in CRPropa 1.4. 
*/

#ifndef _CMBLikeScaling_H_
#define _CMBLikeScaling_H_
#include <TIRBzEvolutionModel.h>

//2. The case of CRPropa 1.4 cmb-like scaling pow(1+z,3) with MaxRedshift at which the IRB will be switched off.
class CMBLikeScaling : public TIRBzEvolutionModel{
 public:
  CMBLikeScaling(std::string XMLFileName);
  double GetScalingFactor(double redshift);
  virtual ~CMBLikeScaling() { }

};
#endif
