#include "utils.h"
#include "crpropaoutput.h"
#include "graphics.h"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>

#include "TObject.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TString.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TList.h"
#include "TCollection.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"
#include "TSystem.h"
#include "TStyle.h"
#include "THStack.h"
#include "TLegend.h"

using namespace std;
/*
void ReadCRs(char* dirname, vector<CRPropaOutput*>& output);
void ReadNeutrinos(char* dirname, vector<CRPropaOutput*>& output);
void ReadPhotons(char* dirname, string& options, vector<CRPropaOutput*>& output);
TH1* ComputeCRSpectrum(const vector<CRPropaOutput*>& output, double spect_ind, double old_spect, int nbinsCR, double EminPr, double Emax, double factor);
TH1* ComputeNeutrinoSpectrum(const vector<CRPropaOutput*>& output, double spect_ind, double old_spect, int nbinsCR, double EminPr, double Emax, double factor, bool& found_neut);
TH1* ComputePhotonSpectrum(const vector<CRPropaOutput*>& output, double spect_ind, double old_spect, int nbins, double Emin, double Emax, double factor, bool& found_phot);
TH1* ComputeMassSpectrum(const vector<CRPropaOutput*>& output, double spect_ind, double old_spect, int nbinsCR, double EminPr, double Emax, double factor);
*/

void plot_spectrum(char* dirname /* Name of single file, or of directory containing several files */, 
		   double norm_en /** Normalization energy, in eV */, 
		   double spect_ind /* New spectral index */, 
		   double old_spect /* Old spectral index */, 
		   double reweightfactor = 1 /* Constant additional reweight factor for pair production in 3D */,
		   double Emax_ = 1e23, // eV; Maximum Energy
		   double EminCR = 1e17, // eV; Minimum Energy for UHECR plot
		   double EminSec = 1e6, // eV; Minimum Energy for gamma-ray and neutrino plots
		   int nbinsCR = 50 /* Number of bins of UHECR plot */, 
		   int nbins = 180 /* Number of bins of gamma-ray plot */
		   ) {

  // Initialization: parameters
  const double Emax = log10(Emax_); // eV; Maximum Energy
  const double EminPr = log10(EminCR); // eV; Minimum Energy for UHECR plot
  const double Emin = log10(EminSec); // eV; Minimum Energy for gamma-ray and neutrino plots

  const double factor    = 1.0; // Choose what to plot. 1 = E^2 dN/dE
  
  InitGraphics();

  // Read AUGER data
  TGraphErrors* grerr = readAugerSpectrum();

  // Read FERMI data
  TGraphAsymmErrors* grFERMI = read_FERMI();

  // Read simulated data found in directory <dirname> : several files can be read at once

  string options = "";
  vector<CRPropaOutput*> output;
  ReadCRs(dirname, output);
  ReadNeutrinos(dirname, output);
  ReadPhotons(dirname, options, output);

  // Filling histos

  bool found_phot = false;
  bool found_neut = false;

  TH1D* prot = (TH1D*)ComputeCRSpectrum(output, "prot", spect_ind, old_spect, nbinsCR, EminPr, Emax, factor);
  TProfile* ms = (TProfile*)ComputeMassSpectrum(output, "ms", spect_ind, old_spect, nbinsCR, EminPr, Emax);
  TH1D* neutrinos = (TH1D*)ComputeNeutrinoSpectrum(output, "neutrinos", spect_ind, old_spect, nbins, Emin, Emax, factor, found_neut);
  TH1D* photons = (TH1D*)ComputePhotonSpectrum(output, "photons", spect_ind, old_spect, nbins, Emin, Emax, factor, reweightfactor, found_phot);

  // Normalizing flux
  const double normEn = log10(norm_en); 
  const double norm = getNorm(grerr, prot, normEn); // Normalize to AUGER data at energy = norm_en

  prot->Scale(norm);
  if (found_phot) photons->Scale(norm);
  if (found_neut) neutrinos->Scale(norm);

  // Drawing
  TCanvas* c2 = new TCanvas("c2", "UHECR spectrum");
  prot->Draw("P");
  grerr->Draw("P");
  gPad->SetLogy();

  if (found_phot) {
    TCanvas* c3 = new TCanvas("c3", "UHE gamma-ray spectrum");
    photons->DrawCopy(options.c_str());
    if (found_phot) grFERMI->Draw("P");
    gPad->SetLogy();
  }

  if (found_neut) {
    TCanvas* c4 = new TCanvas("c4", "UHE neutrino spectrum");
    neutrinos->DrawCopy("P");
    gPad->SetLogy();
  }

  TCanvas* c5 = new TCanvas("c5", "Average mass");
  ms->Draw();

  TCanvas* c6 = new TCanvas("c6", "Combined");
  
  THStack* hst = new THStack("hst", ";log_{10}(E / eV);E^{2}dN/dE (eV m^{-2} s^{-1} sr^{-1})");
  hst->Add(prot);
  if (found_phot) {
    photons->SetMarkerColor(kBlue);
    hst->Add(photons, options.c_str());
  }
  if (found_neut) {
    neutrinos->SetMarkerColor(kGreen);
    hst->Add(neutrinos);
  }
  hst->Draw("NOSTACK");
  grerr->Draw("P");
  if (found_phot) grFERMI->Draw("P");

  gPad->SetLogy();
 
  TLegend* leg1 = new TLegend(0.75, 0.75, 0.95, 0.95);
  leg1->AddEntry(prot, "UHECRs", "lp");
  if (found_neut) leg1->AddEntry(neutrinos, "#nu (all flavors)", "lp");
  if (found_phot) leg1->AddEntry(photons, "$gamma-rays", "lp");
  leg1->AddEntry(grerr, "AUGER data", "lp");
  if (found_phot) leg1->AddEntry(grFERMI, "EBL - FERMI 2010", "lp");

  leg1->Draw();

  for (std::vector<CRPropaOutput*>::iterator i = output.begin(); i != output.end(); ++i) delete *i;
  output.clear();

  return ;
}
