/**
   @class  Kneiske2004_BestFit
   
   @author  Nils Nierstenhoefer nierstenhoefer@physik.uni-wuppertal.de
   
   @brief  Model to scale the mean free path (MFP) for a particle reaction with photons of the infrared background light (IRB) at redshift z=0 to estimate the MFP at a given redshift z'.
   
   The plots from fig.1 of Kneiske, Mannheim, Hartmann (astroph 0309141) show the evolution of the COMOVING MRF spectrum for several reshifts where the wavelength 
   scale corresponds to a COMOVING reference frame. These plots were red from the paper and are converted into 8 TGraphs for z=0, 0.2, 0.4, . . . 4. The units are changed to photon number 
   density [GeV^2] versus photon energy [GeV] (both still in comoving coordinates and h_bar=c=1) and for each of the 8 redshift values the integral is calculated using TGraph::Integral(). 
   These integral values \f$S(z)\f$ are still COMOVING and thus become multiplied with \f$(1+z)^3\f$. After a normalization that guarantees S'(z=0)=1, S(z) versus z is stored
   in a new TGraph and GetScalingFactor(double redshift) returns TGraph::Eval(double ) (using TSpline3). 

   Note, there are only plots available up to z=4. After that value this model returns 0. 

   Cleary the best solution would be take get hold of the Kneiske parametrization to avoid using the graphs. 
*/
#ifndef _Kneiske2004_BestFit_H_
#define _Kneiske2004_BestFit_H_

#include <TIRBzEvolutionModel.h>
#include <CMB.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TAxis.h>

class Kneiske2004_BestFit : public TIRBzEvolutionModel{
 public:
  Kneiske2004_BestFit(std::string XMLFileName);
  double GetScalingFactor(double redshift);
  virtual ~Kneiske2004_BestFit() { if (fIRBScaleFactorTG) delete fIRBScaleFactorTG; }

 protected:
  TGraph* fIRBScaleFactorTG;
};

#endif
