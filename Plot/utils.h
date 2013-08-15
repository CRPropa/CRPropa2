#ifndef _UTILS_H
#define _UTILS_H

#include <iostream>
#include <fstream>
#include <vector>

class TH1D;
class TGraphErrors;
class TGraphAsymmErrors;

double getNorm(TGraphErrors*, TH1D*, const double&);
TGraphErrors* readAugerSpectrum();
TGraphAsymmErrors* read_FERMI();

#endif
