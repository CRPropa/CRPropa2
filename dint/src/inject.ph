#ifndef __inject_ph__
#define __inject_ph__

#ifdef __cplusplus
extern "C" {
#endif

void SetInjectionSpectrum (const PARTICLE part, const dCVector* pEnergy,			  const dCVector* pEnergyWidth, Spectrum* pQ_0,			  const double InjEnergy);

#ifdef __cplusplus
}
#endif

#endif

