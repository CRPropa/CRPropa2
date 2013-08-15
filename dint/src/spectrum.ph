#ifndef __spectrum_ph__
#define __spectrum_ph__

#ifdef __cplusplus
extern "C" {
#endif

void NewSpectrum (Spectrum* pSpectrum, const int num_bins);
void DeleteSpectrum (Spectrum* pSpectrum);
void InitializeSpectrum (Spectrum* pSpectrum);
void InitializeEMSpectrum (Spectrum* pSpectrum);
void InitializeNucleonSpectrum (Spectrum* pSpectrum);
void InitializeNeutrinoSpectrum (Spectrum* pSpectrum);
void SetSpectrum (Spectrum* pSpectrum1, const Spectrum* pSpectrum2);
void SetEMSpectrum (Spectrum* pSpectrum1, const Spectrum* pSpectrum2);
void SetNucleonSpectrum (Spectrum* pSpectrum1, const Spectrum* pSpectrum2);
void SetNeutrinoSpectrum (Spectrum* pSpectrum1, const Spectrum* pSpectrum2);
double GetNumber (const Spectrum* pSpectrum);
double GetEMNumber (const Spectrum* pSpectrum);
double GetNucleonNumber (const Spectrum* pSpectrum);
double GetNeutrinoNumber (const Spectrum* pSpectrum);
double GetEnergy (const Spectrum* pSpectrum, const dCVector* pEnergy);
double GetEMEnergy (const Spectrum* pSpectrum, const dCVector* pEnergy);
double GetNucleonEnergy (const Spectrum* pSpectrum, const dCVector* pEnergy);
double GetNeutrinoEnergy (const Spectrum* pSpectrum, const dCVector* pEnergy);
void DumpSpectrum (const dCVector* pEnergy, const dCVector* pEnergyWidth,		  const Spectrum* pSpectrum, const char* filename);

#ifdef __cplusplus
}
#endif

#endif

