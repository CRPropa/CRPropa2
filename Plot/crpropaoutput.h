#ifndef _CRPROPAOUTPUT_H_
#define _CRPROPAOUTPUT_H_

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>

#include "chealpix.h"
#include "fitsio.h"

#include "TVector3.h"
#include "TLorentzVector.h"
#include "TH1D.h"
#include "TProfile.h"
#include "TObject.h"
#include "TString.h"
#include "TList.h"
#include "TCollection.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"
#include "TSystem.h"

typedef enum {oneD=0, threeD=1} Geometry;


class UHECREvent;
class UHECRNeutrino;
class UHECRPhoton;
class UHECRTrajectory;
class CRPropaOutput;
class CRPropaEvent;
class CRPropaTrajectory;
class CRPropaNeutrino;
class CRPropaPhotonSpectrum;

void ReadCRs(char* dirname, std::vector<CRPropaOutput*>& output);
void ReadNeutrinos(char* dirname, std::vector<CRPropaOutput*>& output);
void ReadPhotons(char* dirname, std::string& options, std::vector<CRPropaOutput*>& output);
TH1* ComputeCRSpectrum(const std::vector<CRPropaOutput*>& output, char* histoname, double spect_ind, double old_spect, int nbinsCR, double EminPr, double Emax, double factor);
TH1* ComputeMassSpectrum(const std::vector<CRPropaOutput*>& output, char* histoname, double spect_ind, double old_spect, int nbinsCR, double EminPr, double Emax);
TH1* ComputeNeutrinoSpectrum(const std::vector<CRPropaOutput*>& output, char* histoname, double spect_ind, double old_spect, int nbinsCR, double EminPr, double Emax, double factor, bool& found_neut);
TH1* ComputePhotonSpectrum(const std::vector<CRPropaOutput*>& output, char* histoname, double spect_ind, double old_spect, int nbins, double Emin, double Emax, double factor, double reweightfactor, bool& found_phot);
TProfile* ComputeCumulativeDeflections(const std::vector<CRPropaOutput*>& output, const char* histoname, const double& spect_ind, const double& old_spect, const int& nbins, const double& Emin, const double& Emax, const double& Dmax);
TProfile* ComputeDeflections(const std::vector<CRPropaOutput*>& output, const char* histoname, const double& spect_ind, const double& old_spect, const int& nbins, const double& Emin, const double& Emax);
float* ComputeSkyMap(const std::vector<CRPropaOutput*>& output, const double& newindex, const double& oldindex, const double& Emin, const double& Emax, const double& Dmax, const long& nside);

class UHECREvent {
 protected:
  int type;
  int inittype;
  TVector3 inpos;
  TLorentzVector inmom;
  double time;
  TVector3 pos;
  TLorentzVector mom;

 public:
  UHECREvent() { }
  UHECREvent(int type_, int inittype_, TVector3 inpos_, TLorentzVector inmom_, double time_, TVector3 pos_, TLorentzVector mom_) :
    type(type_),
    inittype(inittype_),
    inpos(inpos_),
    inmom(inmom_),
    time(time_),
    pos(pos_),
    mom(mom_)
      { }

  ~UHECREvent() { }

  inline int GetType() { return type; }
  inline int GetInType() { return inittype; } 
  inline TVector3 GetInPos() { return inpos; }
  inline TLorentzVector GetInMom() { return inmom; }
  inline double GetTime() { return time; }
  inline TVector3 GetPos() { return pos; }
  inline TLorentzVector GetMom() { return mom; }
  inline double GetInEnergy() { return inmom.T(); }
  inline double GetEnergy() { return mom.T(); }
};

class UHECRNeutrino {
 protected:
  int type;
  int parenttype;
  TVector3 sourcepos;
  TVector3 inpos;
  double inE;
  TVector3 pos;
  TLorentzVector mom;
  double sourceE;

 public:
  UHECRNeutrino() { }
  UHECRNeutrino(int type_, int parenttype_, TVector3 sourcepos_, TVector3 inpos_, double inE_, TVector3 pos_, TLorentzVector mom_, double sourceE_) :
    type(type_),
    parenttype(parenttype_),
    sourcepos(sourcepos_),
    inpos(inpos_),
    inE(inE_),
    pos(pos_),
    mom(mom_),
    sourceE(sourceE_)
      { }

  ~UHECRNeutrino() { }

  inline int GetType() { return type; }
  inline int GetParentType() { return parenttype; } 
  inline TVector3 GetInPos() { return inpos; }
  inline TVector3 GetSourcePos() { return sourcepos; }
  inline TVector3 GetPos() { return pos; }
  inline TLorentzVector GetMom() { return mom; }
  inline double GetInEnergy() { return inE; }
  inline double GetSourceEnergy() { return sourceE; }
  inline double GetEnergy() { return mom.T(); }
};

class UHECRPhoton {
 protected:
  int type;
  int inittype;
  std::string cascorig;
  TVector3 sourcepos;
  TVector3 inpos;
  TVector3 pos;
  TLorentzVector mom;
  double sourceE;
  std::vector<double> spectrum;

 public:
  UHECRPhoton() { }
  UHECRPhoton(int type_, int inittype_, std::string cascorig_, TVector3 sourcepos_, TVector3 inpos_, TVector3 pos_, TLorentzVector mom_, double sourceE_, std::vector<double> spectrum_) :
    type(type_),
    inittype(inittype_),
    cascorig(cascorig_),
    sourcepos(sourcepos_),
    inpos(inpos_),
    pos(pos_),
    mom(mom_),
    sourceE(sourceE_),
    spectrum(spectrum_)
      { }

  ~UHECRPhoton() { }

  inline int GetType() { return type; }
  inline int GetInType() { return inittype; } 
  inline std::string GetCascadeOrig() { return cascorig; }
  inline TVector3 GetInPos() { return inpos; }
  inline TVector3 GetSourcePos() { return sourcepos; }
  inline TVector3 GetPos() { return pos; }
  inline TLorentzVector GetMom() { return mom; }
  inline double GetSourceEnergy() { return sourceE; }
  inline double GetEnergy() { return mom.T(); }
  inline std::vector<double> GetSpectrum() { return spectrum; }

};

class UHECRTrajectory {
 protected:
  int type;
  int inittype;
  double time;
  TVector3 pos;
  double En;

 public:
  UHECRTrajectory() { }
  UHECRTrajectory(int type_, int inittype_, double time_, TVector3 pos_, double En_) :
    type(type_),
    inittype(inittype_),
    time(time_),
    pos(pos_),
    En(En_)
      { }

  ~UHECRTrajectory() { }

  inline int GetType() { return type; }
  inline int GetInType() { return inittype; } 
  inline double GetTime() { return time; }
  inline TVector3 GetPos() { return pos; }
  inline double GetEnergy() { return En; }
};


class CRPropaOutput : public TObject {
  
 protected:
  Geometry geo;
  std::string whoami;
  bool Ok;

 public:
  virtual ~CRPropaOutput() {}

  inline std::string WhoAmI() { return whoami; }
  inline Geometry GetGeo() { return geo; }
  inline bool IsOk() { return Ok; }

  virtual TH1D* ComputeSpectrum(const char*, const double&, const double&, const int&, const double&, const double&, const double&, const double& reweightfactor = 1.0) { return NULL; }
  virtual TProfile* ComputeMassSpectrum(const char*, const double&, const double&, const int&, const double&, const double&) { return NULL; }
  virtual TProfile* ComputeDeflections(const char*, const double&, const double&, const int&, const double&, const double&) { return NULL; }
  virtual TProfile* ComputeCumulativeDeflections(const char*, const double&, const double&, const int&, const double&, const double&, const double&) { return NULL; }
  virtual std::vector<double> ComputeSkyMap(const double&, const double&, const double&, const double&, const double&, const long&) { return std::vector<double>(); }
  virtual void SkyMapEvents(const char*, const double&, const double&, const double&) { /*...*/ }
  virtual void SkyMapSources(const char*, const double&) { /*...*/ }
  virtual void DrawTrajectory(unsigned int) { /*...*/ }
  
  ClassDef(CRPropaOutput, 1);
};

class CRPropaEvent : public CRPropaOutput {

 public :
  CRPropaEvent() {}
  CRPropaEvent(std::string&, const std::string&);
  CRPropaEvent(const Geometry& geo_) { geo = geo_; }

  virtual ~CRPropaEvent() { for (std::vector<UHECREvent*>::iterator i = ev.begin(); i != ev.end(); ++i) delete *i; ev.clear(); }

  inline std::vector<UHECREvent*> GetArray() { return ev; }

  void Add(CRPropaEvent* cr) {
    if (geo != cr->GetGeo()) { Error("Add", "Trying to add UHECRs inconsistently. geo = %i, cr->GetGeo() = %i", geo, cr->GetGeo()); return ; }

    std::vector<UHECREvent*> ev1 = cr->GetArray();
    for (size_t i = 0; i < ev1.size(); ++i) ev.push_back(ev1[i]);
    return ;
  }
  

  virtual TH1D* ComputeSpectrum(const char*, const double&, const double&, const int&, const double&, const double&, const double&, const double& reweightfactor = 1.0);
  virtual TProfile* ComputeMassSpectrum(const char*, const double&, const double&, const int&, const double&, const double&);
  virtual TProfile* ComputeDeflections(const char*, const double&, const double&, const int&, const double&, const double&);
  virtual TProfile* ComputeCumulativeDeflections(const char*, const double&, const double&, const int&, const double&, const double&, const double&);
  virtual std::vector<double> ComputeSkyMap(const double&, const double&, const double&, const double&, const double&, const long&);
  virtual void SkyMapEvents(const char*, const double& Emin = 0, const double& Emax = 0, const double& Dmax = 0);
  virtual void SkyMapSources(const char*, const double& Dmax = 0);

 protected :
  std::vector<UHECREvent*> ev;

  ClassDef(CRPropaEvent, 1);
};

class CRPropaTrajectory : public CRPropaOutput {

 public :
  CRPropaTrajectory() {}
  CRPropaTrajectory(std::string&, const std::string&);

  virtual ~CRPropaTrajectory() { for (std::vector<UHECRTrajectory*>::iterator i = traj.begin(); i != traj.end(); ++i) delete *i; traj.clear(); ntracks.clear(); }
  virtual void DrawTrajectory(unsigned int);

 protected :
  std::vector<UHECRTrajectory*> traj;
  std::vector<int> ntracks;
    

  ClassDef(CRPropaTrajectory, 1);
};


class CRPropaNeutrino : public CRPropaOutput {

 public :
  CRPropaNeutrino() {}
  CRPropaNeutrino(std::string&, const std::string&);

  virtual ~CRPropaNeutrino() { for (std::vector<UHECRNeutrino*>::iterator i = neu.begin(); i != neu.end(); ++i) delete *i; neu.clear(); }
  virtual TH1D* ComputeSpectrum(const char*, const double&, const double&, const int&, const double&, const double&, const double&, const double& reweightfactor = 1.0);

 protected :
  std::vector<UHECRNeutrino*> neu;

  ClassDef(CRPropaNeutrino, 1);
};

class CRPropaPhotonSpectrum : public CRPropaOutput {

 public :
  CRPropaPhotonSpectrum() {}
  CRPropaPhotonSpectrum(std::string&, const std::string&);

  virtual ~CRPropaPhotonSpectrum() { for (std::vector<UHECRPhoton*>::iterator i = phot.begin(); i != phot.end(); ++i) delete *i; phot.clear(); }
  virtual TH1D* ComputeSpectrum(const char*, const double&, const double&, const int&, const double&, const double&, const double&, const double& reweightfactor = 1.0);
  virtual void SkyMapEvents(const char*, const double&, const double&, const double&) { /*...*/ }
  virtual void SkyMapSources(const char*, const double& Dmax = 0.0);

 protected :
  std::vector<UHECRPhoton*> phot;
  static const double EnergyBins[170];

  ClassDef(CRPropaPhotonSpectrum, 1);
};

const double CRPropaPhotonSpectrum::EnergyBins[170] = { 1.0066347e-11, 1.2672780e-11, 1.5954085e-11, 2.0085003e-11,
							2.5285520e-11, 3.1832584e-11, 4.0074849e-11, 5.0451246e-11,
							6.3514355e-11, 7.9959836e-11, 1.0066347e-10, 1.2672780e-10,
							1.5954085e-10, 2.0085003e-10, 2.5285520e-10, 3.1832584e-10,
							4.0074849e-10, 5.0451246e-10, 6.3514355e-10, 7.9959836e-10,
							1.0066347e-09, 1.2672780e-09, 1.5954085e-09, 2.0085003e-09,
							2.5285520e-09, 3.1832584e-09, 4.0074849e-09, 5.0451246e-09,
							6.3514355e-09, 7.9959836e-09, 1.0066347e-08, 1.2672780e-08,
							1.5954085e-08, 2.0085003e-08, 2.5285520e-08, 3.1832584e-08,
							4.0074849e-08, 5.0451246e-08, 6.3514355e-08, 7.9959836e-08,
							1.0066347e-07, 1.2672780e-07, 1.5954085e-07, 2.0085003e-07,
							2.5285520e-07, 3.1832584e-07, 4.0074849e-07, 5.0451246e-07,
							6.3514355e-07, 7.9959836e-07, 1.0066347e-06, 1.2672780e-06,
							1.5954085e-06, 2.0085003e-06, 2.5285520e-06, 3.1832584e-06,
							4.0074849e-06, 5.0451246e-06, 6.3514355e-06, 7.9959836e-06,
							1.0066347e-05, 1.2672780e-05, 1.5954085e-05, 2.0085003e-05,
							2.5285520e-05, 3.1832584e-05, 4.0074849e-05, 5.0451246e-05,
							6.3514355e-05, 7.9959836e-05, 0.00010066347, 0.00012672780,
							0.00015954085, 0.00020085003, 0.00025285520, 0.00031832584,
							0.00040074849, 0.00050451246, 0.00063514355, 0.00079959836,
							0.0010066347,  0.0012672780,  0.0015954085,  0.0020085003,
							0.0025285520,  0.0031832584,  0.0040074849,  0.0050451246,
							0.0063514355,  0.0079959836,  0.010066347,   0.012672780,
							0.015954085,   0.020085003,   0.025285520,   0.031832584,
							0.040074849,   0.050451246,   0.063514355,   0.079959836,
							0.10066347,    0.12672780,    0.15954085,    0.20085003,
							0.25285520,    0.31832584,    0.40074849,    0.50451246,
							0.63514355,    0.79959836,    1.0066347,     1.2672780,
							1.5954085,     2.0085003,     2.5285520,     3.1832584,
							4.0074849,     5.0451246,     6.3514355,     7.9959836,
							10.066347,     12.672780,     15.954085,     20.085003,
							25.285520,     31.832584,     40.074849,     50.451246,
							63.514355,     79.959836,     100.66347,     126.72780,
							159.54085,     200.85003,     252.85520,     318.32584,
							400.74849,     504.51246,     635.14355,     799.59836,
							1006.6347,     1267.2780,     1595.4085,     2008.5003,
							2528.5520,     3183.2584,     4007.4849,     5045.1246,
							6351.4355,     7995.9836,     10066.347,     12672.780,
							15954.085,     20085.003,     25285.520,     31832.584,
							40074.849,     50451.246,     63514.355,     79959.836,
							100663.47,     126727.80,     159540.85,     200850.03,
							252855.20,     318325.84,     400748.49,     504512.46,
							635143.55,     799598.36}; // in EeV
#endif
