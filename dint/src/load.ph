#ifndef __load_ph__
#define __load_ph__

#ifdef __cplusplus
extern "C" {
#endif

void LoadICSTables (RawTotalRate* ICSTotalRate, RawDiffRate* ICSPhotonRate,		   RawDiffRate* ICSScatRate, const int num_main_bins,		   char* aDirTables);
void LoadPPTables (RawTotalRate* PPTotalRate, RawDiffRate* PPDiffRate, 		  const int num_main_bins, char* aDirTables);
void LoadTPPTables (RawTotalRate* TPPTotalRate, RawDiffRate* TPPDiffRate, 		   const int num_main_bins, char* aDirTables);
void LoadDPPTables (RawTotalRate* DPPRate, const int num_main_bins,		   char* aDirTables);
void LoadPPPNucleonTables (RawTotalRate* PPPProtonLossRate, 			  RawTotalRate* PPPNeutronLossRate,			  RawDiffRate* PPPProtonScatRate,			  RawDiffRate* PPPProtonNeutronRate,			  RawDiffRate* PPPNeutronProtonRate,			  const int num_main_bins,			  char* aDirTables);
void LoadPPPEMTables (RawDiffRate* PPPProtonPhotonRate,		     RawDiffRate* PPPProtonElectronRate,		     RawDiffRate* PPPProtonPositronRate,		     RawDiffRate* PPPNeutronElectronRate, 		     const int num_main_bins,		     char* aDirTables);
void LoadNPPNucleonTables (RawTotalRate* NPPTotalRate, const int num_main_bins,			  char* aDirTables);
void LoadNPPSecondaryTables (RawDiffRate* NPPDiffRate, const int num_main_bins,			    char* aDirTables);
void LoadPPPNeutrinoTables (RawDiffRate* PPPProtonElectronNeutrinoRate,			   RawDiffRate* PPPProtonAntiElectronNeutrinoRate,			   RawDiffRate* PPPProtonMuonNeutrinoRate,			   RawDiffRate* PPPProtonAntiMuonNeutrinoRate,			   RawDiffRate* PPPNeutronAntiElectronNeutrinoRate,			   RawDiffRate* PPPNeutronMuonNeutrinoRate,			   RawDiffRate* PPPNeutronAntiMuonNeutrinoRate, 			   const int num_main_bins, char* aDirTables);
void LoadNeutronDecayNucleonTables (TotalRate* neutronDecayRate,				   DiffRate* neutronDecayProtonRate,				   const int num_main_bins,				   char* aDirTables);
void LoadNeutronDecaySecondaryTables (DiffRate* neutronDecayElectronRate, 				     const int num_main_bins,				     char* aDirTables);
void LoadNeutrinoTables (const int tauNeutrinoMassSwitch,			TotalRate* NNElNeutTotalRate, 			TotalRate* NNMuonNeutTotalRate,			TotalRate* NNTauNeutTotalRate, 			DiffRate* NNElNeutScatRate, 			DiffRate* NNElNeutMuonNeutRate,			DiffRate* NNElNeutTauNeutRate, 			DiffRate* NNElNeutElectronRate, 			DiffRate* NNElNeutPhotonRate,			DiffRate* NNElNeutProtonRate, 			DiffRate* NNMuonNeutScatRate, 			DiffRate* NNMuonNeutElNeutRate,			DiffRate* NNMuonNeutTauNeutRate, 			DiffRate* NNMuonNeutElectronRate, 			DiffRate* NNMuonNeutPhotonRate, 			DiffRate* NNMuonNeutProtonRate, 			DiffRate* NNTauNeutScatRate,			DiffRate* NNTauNeutElNeutRate, 			DiffRate* NNTauNeutMuonNeutRate, 			DiffRate* NNTauNeutElectronRate, 			DiffRate* NNTauNeutPhotonRate, 			DiffRate* NNTauNeutProtonRate, 			const int num_main_bins,			char* aDirTables);

#ifdef __cplusplus
}
#endif

#endif

