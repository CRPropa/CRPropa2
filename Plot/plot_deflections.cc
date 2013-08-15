#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>

#include "chealpix.h"
#include "fitsio.h"

#include "TH1D.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TList.h"
#include "TCollection.h"
#include "TObject.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TString.h"

#include "crpropaoutput.h"
#include "utils.h"
#include "graphics.h"

using namespace std;

//vector<CRPropaOutput*> ReadSimulations(char* dirname, string& options);

void plot_deflections(char* dirname, /* Name of single file, or of a directory containing several files */
		      double norm_en /* Normalization energy, in eV */, 
		      char* mapn /* Name of Skymap file */,
		      double spect_ind, /* New spectral index */
		      double old_spect, /* Old spectral index */
		      double Dmax, /* Mpc; maximum distance of incoming UHECR */
		      double Emax_ = 1e23, /* eV; Maximum energy of UHECR plot */
		      double EminCR = 1e17,  /* eV; Minimum energy of UHECR plot */
		      int nbinsPr = 100, /* Number of bins of UHECR plot */
		      int nbinsDefl = 90, /* Number of bins of deflection plot */
		      double EminDefl = 5,  /* EeV; minimal energy to compute deflections */
		      double EmaxDefl = 200, /* EeV; maximal energy to compute deflections */
		      double EminSM = 60,  /* EeV; minimal energy to compute sky map */
		      double EmaxSM = 200 /* EeV; maximal energy to compute sky map */
		      ) {
  
  TString checkmapn(mapn);
  if (!checkmapn.Contains(".fits")) {
    cerr << "Skymap should be in FITS format" << endl;
    return ;
  }

  // Initialization: parameters
  const double Emax = log10(Emax_); // eV; Maximum energy of UHECR plot
  const double EminPr = log10(EminCR); // eV; Minimum energy of UHECR plot

  const double factor = 1.0; // Overall normalization factor

  // Sky map parameters, see Healpix
  const int resolution = 4;
  const long nside = (long)pow(2.0, resolution);
  const long npix = nside2npix(nside);

  // Initialization: graphs and plots
  InitGraphics();

  // Read AUGER data
  TGraphErrors* grerr = readAugerSpectrum();

  // Read simulated data found in directory <dirname> : several files can be read at once

  string options = "";
  vector<CRPropaOutput*> output;
  ReadCRs(dirname, output);
  
  // Filling histos
  //  TH1D* protLI = new TH1D("protLI", ";log_{10}(E / eV);E^{2}dN/dE (eV m^{-2} s^{-1} sr^{-1})", nbinsPr, EminPr, Emax); 
  //TProfile* cumdefl = new TProfile("cumdefl", ";#delta (#circ);", nbinsDefl, 0, 180); 
  //TProfile* deflections = new TProfile("deflections", ";log_{10}(E / eV);#delta (#circ)", nbinsPr, EminPr, Emax); 

  //string name;

  //TH1D* h4;
  //TProfile* h3;
  TH1D* prot = (TH1D*)ComputeCRSpectrum(output, "prot", spect_ind, old_spect, nbinsPr, EminPr, Emax, factor);
  TProfile* cumdefl = (TProfile*)ComputeCumulativeDeflections(output, "cumdefl", spect_ind, old_spect, nbinsDefl, EminDefl, EmaxDefl, Dmax);;
  TProfile* deflections = (TProfile*)ComputeDeflections(output, "deflections", spect_ind, old_spect, nbinsPr, EminPr, Emax);
  float* skymap = ComputeSkyMap(output, spect_ind, old_spect, EminSM, EmaxSM, Dmax, nside);

  /*
  for (unsigned int i = 0; i < output.size(); i++) {
    if (output[i]->IsOk()) {
      if (output[i]->WhoAmI() == "CR") {
	name ="h4";
	name += i;
	h4 = output[i]->ComputeSpectrum(name.c_str(), spect_ind, old_spect, nbinsPr, EminPr, Emax, factor);
	protLI->Add(h4);
	h4->Delete();
	
	name ="h3";
	name += i;
	h3 = output[i]->ComputeCumulativeDeflections(name.c_str(), spect_ind, old_spect, nbinsDefl, EminDefl, EmaxDefl, Dmax);
	cumdefl->Add(h3);
	h3->Delete();
	
	vector<double> fakemap = output[i]->ComputeSkyMap(spect_ind, old_spect, EminSM, EmaxSM, Dmax, nside);
	for (long j = 0; j < npix; ++j) {
	  skymap[j] += fakemap[j];
	}
	
	name = "h6";
	name += i;
	h3 = output[i]->ComputeDeflections(name.c_str(), spect_ind, old_spect, nbinsPr, EminPr, Emax);
	deflections->Add(h3);
	h3->Delete();
      }
    }
  }
  */

  // Normalizing skymap
  float* maximum = max_element(skymap, skymap+npix);
  for (int i = 0; i < npix; ++i) skymap[i] /= (*maximum);

  string mapname("!");
  mapname += mapn;

  // Plot skymap in Galactic coordinates
  if (write_healpix_map(skymap, nside, (char*)mapname.c_str(), '0', "G")) cerr << "Something wrong with writing the skymap." << endl; 

  // Normalizing flux
  const double normEn = log10(norm_en);
  const double norm = getNorm(grerr, prot, normEn);

  prot->Scale(norm);

  // Drawing
  TCanvas* c2 = new TCanvas("c2", " ");
  prot->Draw();
  grerr->Draw("P");
  gPad->SetLogy();


  TCanvas* c3 = new TCanvas("c3", "Deflections");
  deflections->Draw();

  TCanvas* c4 = new TCanvas("c4", "Cumulative deflections");
  cumdefl->Draw();

  for (std::vector<CRPropaOutput*>::iterator i = output.begin(); i != output.end(); ++i) delete *i;
  output.clear();

  return ;
}

