#ifndef __sync_ph__
#define __sync_ph__

#ifdef __cplusplus
extern "C" {
#endif

void LoadSyncTable (dCVector* syncTable);
void InitializeSynchrotron (const double B_0, const dCVector* pEnergy, 			   const dCVector* pEnergyWidth,			   dCVector* synchrotronLoss, DiffRate* syncRate);

#ifdef __cplusplus
}
#endif

#endif

