#ifndef __background_ph__
#define __background_ph__

#ifdef __cplusplus
extern "C" {
#endif

void LoadPhotonBackground (const double redshift,                           dCVector* pBgEnergy, dCVector* pBgEnergyWidth,                           dCVector* pBgPhotonDensity);
void LoadCMB (const double redshift, const dCVector* pBgEnergy,              const dCVector* pBgEnergyWidth, dCVector* pBgPhotonDensity);
void LoadIR (const double redshift, const dCVector* pBgEnergy,            const dCVector* pBgEnergyWidth, dCVector* pBgPhotonDensity);
double HighIR (const double zTarget, const double zObserve, 	      const double energy0, const double deltaO,	      const double deltaD);
double LowIR (const double zTarget, const double zObserve, 	     const double energy0, const double deltaO,	     const double deltaD);
double OpticalIR (const double energy);
double DustIR (const double energy);
void LoadRadio (const double redshift, const dCVector* pBgEnergy,               const dCVector* pBgEnergyWidth, dCVector* pBgPhotonDensity);
double HighRadio (const double zTarget, const double zObserve, 		const double energy0);
double MedRadio (const double zTarget, const double zObserve, 		const double energy0);
double ObsRadio (const double zTarget, const double zObserve, 		const double energy0);

#ifdef __cplusplus
}
#endif

#endif

