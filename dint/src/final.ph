#ifndef __final_ph__
#define __final_ph__

#ifdef __cplusplus
extern "C" {
#endif

void CheckEnergy (const int sourceTypeSwitch, const double brightPhaseExp,                 const double startingRedshift, 		 const double rightRedshift, const Spectrum* pSpectrum, 		 const dCVector* pEnergy, const double initialTotalEnergy);
void FinalPrintOutToTheScreen (const double distance, 			      const double startingRedshift,			      const double propagatingDistance);

#ifdef __cplusplus
}
#endif

#endif

