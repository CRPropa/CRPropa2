#ifndef __decay_ph__
#define __decay_ph__

#ifdef __cplusplus
extern "C" {
#endif

double PionToPhoton (const int iPhoton, const int iPion);
double PionToLepton (const double leptonEnergy, const double pionEnergy);
double PionToElectronNeutrino (const double neutrinoEnergy, 			      const double pionEnergy);
double PionToMuonNeutrino (const int iNeutrino, const int iPion);

#ifdef __cplusplus
}
#endif

#endif

