#include "utils.h"

#define BAD_NORM_ENERGY -2

#include "TGraphErrors.h"
#include "TH1D.h"
#include "TSpline.h"
#include "TGraphAsymmErrors.h"

using namespace std;

double getNorm(TGraphErrors* grerr, TH1D* prot, const double& normEn) {
  int numbin = prot->FindBin(normEn);
  if (prot->GetBinContent(numbin) == 0) {
    cerr << "You want normalize the UHECR spectrum at E = " << pow(10, normEn) << " eV. However, the UHECR flux from the MonteCarlo is 0 at this energy." << endl;
    exit(BAD_NORM_ENERGY);
  }
  TSpline5 spl("spl", grerr);
  return spl.Eval(normEn)/prot->GetBinContent(numbin);
}

TGraphErrors* readAugerSpectrum() {
  vector<double> energyAuger;
  vector<double> fluxAuger;
  vector<double> errflAuger;
	
  //  ifstream infile("DATA/Auger_Spectrum_2007_data.dat", ios::out);
  ifstream infile("DATA/combined_spectrum_icrc09.txt", ios::out);
  if (!infile.is_open()) {
    cerr << "AUGER spectrum FILE not found!" << endl;
    return NULL;
  }
	
  double var1 = 0;
  double var2 = 0;
  double var3 = 0;
  double var4 = 0;
	
  char s[3000];
  for (int i = 0; i < 22; ++i) infile.getline(s,3000);

   while (infile.good()) {
    infile >> var1 >> var2 >> var3 >> var4;
    energyAuger.push_back(var1);
    fluxAuger.push_back(var2/1e30*pow(10, var1*2));
    errflAuger.push_back(var3/var2*fluxAuger.back());
  }
	
  const int NAuger = energyAuger.size()-1;
	
  vector<double> errenAuger(NAuger,0.0);
	
  double enAug[NAuger];
  double errenAug[NAuger];
  double fluxAug[NAuger];
  double errflAug[NAuger];
	
  for (int i = 0; i < NAuger; i++) {
    enAug[i] = energyAuger[i];
    errenAug[i] = errenAuger[i];
    fluxAug[i] = fluxAuger[i];
    errflAug[i] = errflAuger[i];
  }
  TGraphErrors* gr = new TGraphErrors(NAuger, enAug, fluxAug, errenAug, errflAug);
  gr->SetLineColor(kRed);
  gr->SetMarkerColor(kRed);
  gr->SetMarkerStyle(7);

  return gr;
}

TGraphAsymmErrors* read_FERMI() {

  ifstream infile("DATA/gEBL_Fermi.dat", ios::in);

  vector<double> Ecenter;
  vector<double> ErrElow;
  vector<double> ErrEhigh;
  vector<double> Flux;
  vector<double> ErrFluxlow;
  vector<double> ErrFluxhigh;

  double Emin, Emax, flux, errflux;
  char s[256];

  infile.getline(s, 256);
  while (infile >> Emin >> Emax >> flux >> errflux) {
    Ecenter.push_back(log10(0.5*(Emin+Emax))+9.0);
    ErrEhigh.push_back(9.0+log10(Emax)-Ecenter.back());
    ErrElow.push_back(Ecenter.back()-log10(Emin)-9.0);
    Flux.push_back(1.e4*pow(10.0, Ecenter.back())*flux);
    ErrFluxlow.push_back(1.e4*pow(10.0, Ecenter.back())*errflux);
    ErrFluxhigh.push_back(1.e4*pow(10.0, Ecenter.back())*errflux);
  }
  infile.close();

  

  TGraphAsymmErrors* gr = new TGraphAsymmErrors(Flux.size(), &(Ecenter[0]), &(Flux[0]), &(ErrElow[0]), &(ErrEhigh[0]), &(ErrFluxlow[0]), &(ErrFluxhigh[0]));
  gr->SetMarkerStyle(7);
  gr->SetMarkerColor(kMagenta);
  gr->SetLineColor(kMagenta);

  return gr;

}
