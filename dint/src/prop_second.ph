#ifndef __prop_second_ph__
#define __prop_second_ph__

#ifdef __cplusplus
extern "C" {
#endif

void DumpEnergy (const Spectrum* pSpectrum, const dCVector* pEnergy);
void secondary_propagation (const double dist_observer, const double InjEnergy,			   PARTICLE part, dCVector* pB_field,			   const dCVector* pEnergy,			   const dCVector* pEnergyWidth,			   Spectrum* pSpectrum,			   char* aDirTables);

#ifdef __cplusplus
}
#endif

#endif

