#ifndef __prepare_ph__
#define __prepare_ph__

#ifdef __cplusplus
extern "C" {
#endif

void SetEnergyBins (const int min_energy_exp, dCVector* pEnergy, 		   dCVector* pEnergyWidth);
void SetDeltaG (const dCVector* pEnergy, dCVector* pDeltaG);
void GetModeOfInput (FILE* input, int* pInputMethodSwitch);
void BasicParameterInput (FILE* input, const int argc, int* pMinEnergyExp,			 int* pNumSmallSteps, 			 double* pConvergeParameter,			 double* pStartingRedshift, 			 double* pStartingDistanceInMpc);
void InteractionParameterInput (FILE* input, const int argc,			       int* pSynchrotronSwitch, double* pB_0, 			       int* pTauNeutrinoMassSwitch, int* pICSSwitch, 			       int* pPPSwitch, int* pTPPSwitch, 			       int* pDPPSwitch, int* pPPPSwitch, 			       int* pNPPSwitch, int* pNeutronDecaySwitch,			       int* pNucleonToSecondarySwitch,			       int* pNeutrinoNeutrinoSwitch);
void ModelParameterInput (FILE* input, const int argc,			 int* pSourceTypeSwitch, double* pMinDistance,			 double* pBrightPhaseExp, int* pModelTypeSwitch);
void PrepareSpectra (const int sourceTypeSwitch, const Spectrum* pQ_0,                     Spectrum* pSpectrum, Spectrum* pSpectrumNew, 		    Spectrum* pDerivative);
void ComputeTotalInitialContent (const dCVector* pEnergy, const Spectrum* pQ_0,                                 double* initialPhotonEnergy,                                double* initialLeptonEnergy, 				double* initialNucleonEnergy,				double* initialNeutrinoEnergy,                                double* initialTotalEnergy,                                double* initialPhotonNumber,                                double* initialLeptonNumber,				double* initialNucleonNumber,				double* initialNeutrinoNumber,                                double* initialTotalNumber);
void ComputeContinuousEnergyLoss (const int synchrotronSwitch,                                  const dCVector* synchrotronLoss,                                  const dCVector* otherLoss,                                  dCVector* continuousLoss);

#ifdef __cplusplus
}
#endif

#endif

