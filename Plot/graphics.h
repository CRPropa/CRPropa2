#include "TStyle.h"
#include "TH1.h"

void InitGraphics(void) {
  // Initialization: graphs and plots
  TH1::SetDefaultSumw2(kTRUE);
  gStyle->SetMarkerStyle(8);
  gStyle->SetMarkerSize(0.5);
  gStyle->SetPalette(1);
  gStyle->SetOptStat(000000000);
}
