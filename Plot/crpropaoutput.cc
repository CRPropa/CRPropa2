#include "crpropaoutput.h"

#include "TFile.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TH2D.h"
#include "TMath.h"
#include "TGeoTrack.h"
#include "TParticle.h"
#include "TGeoManager.h"

using namespace std;

/*
*********************************************************************
*********************************************************************
*/

CRPropaEvent::CRPropaEvent(string& filename, const string& filetype) {

  float type = 0;
  float inittype = 0;
  float initx = 0;
  float inity = 0;
  float initz = 0;
  float x = 0;
  float y = 0;
  float z = 0;
  float inite = 0;
  float inittheta = 0;
  float initphi = 0;
  float energy = 0;
  float theta = 0;
  float phi = 0;
  float t = 0;


  Ok = false;

  if (filename == "." || filename == "..") return ;

  if (filetype == "ROOT") {

    TFile rootfile(filename.c_str(), "READ");

    if (rootfile.IsOpen()) {
      TNtuple* t1;
      rootfile.GetObject("events",t1);

      if (t1) {

	Ok = true;
	TString title(t1->GetTitle());
	if (title.Contains("1D")) geo = oneD;
	else geo = threeD;

	t1->SetBranchAddress("Particle_Type", &type);
	t1->SetBranchAddress("Initial_Type", &inittype);
	t1->SetBranchAddress("Time_Mpc", &t);

	if (geo == oneD) {
	  t1->SetBranchAddress("Initial_Position_Mpc", &initx);
	  t1->SetBranchAddress("Energy_EeV", &energy);
	  t1->SetBranchAddress("Initial_Energy_EeV", &inite);
	}
	else {
	  t1->SetBranchAddress("Initial_Position_X_Mpc", &initx);
	  t1->SetBranchAddress("Initial_Position_Y_Mpc", &inity);
	  t1->SetBranchAddress("Initial_Position_Z_Mpc", &initz);

	  t1->SetBranchAddress("Position_X_Mpc", &x);
	  t1->SetBranchAddress("Position_Y_Mpc", &y);
	  t1->SetBranchAddress("Position_Z_Mpc", &z);

	  t1->SetBranchAddress("Initial_Momentum_E_EeV", &inite);
	  t1->SetBranchAddress("Initial_Momentum_theta", &inittheta);
	  t1->SetBranchAddress("Initial_Momentum_phi", &initphi);

	  t1->SetBranchAddress("Momentum_E_EeV", &energy);
	  t1->SetBranchAddress("Momentum_theta", &theta);
	  t1->SetBranchAddress("Momentum_phi", &phi);
	}

	int Nev = t1->GetEntries();

	for (int i = 0; i < Nev; i++) {
	  t1->GetEntry(i);

	  TLorentzVector v1(inite,0,0,inite);
	  v1.SetTheta(inittheta);
	  v1.SetPhi(initphi);

	  TLorentzVector v2(energy,0,0,energy);
	  v2.SetTheta(theta);
	  v2.SetPhi(phi);

	  ev.push_back( new UHECREvent((int)type, (int)inittype, TVector3(initx,inity,initz), v1, t, TVector3(x,y,z), v2) );

	}
	rootfile.Close();
	t1 = NULL;
      } // t1
      else Warning("CRPropaEvent", "Your ROOT file\"%s\" did not contain any UHECR.", filename.c_str());
    } // rootfile
    else Error("CRPropaEvent", "Could not open the ROOT file \"%s\" in CRPropaEvent.", filename.c_str());
  } // strcmp
  else if (filetype == "ASCII") {
    ifstream infile(filename.c_str(), ios::in);
    if (!infile.is_open()) Error("CRPropaEvent", "Could not open your ASCII file \"%s\".", filename.c_str());

    string str("");
    char* s = new char[3000];
    
    infile.getline(s,3000);
    if (!strcmp(s, "#CRPropa - Output data file")) {
      infile.getline(s,3000);

      TString rootstring(s);
      if (rootstring.Contains("[")) geo = threeD;
      else geo = oneD;

      while (infile.good()) {
	infile >> str;

	if (str[0] == '#') {
	  if (str[1] == 'O') break;
	  infile.ignore(256,'\n');
	}
	else {
	  type = atof(str.c_str());
	  if (type > 1000) {
	    Ok = true;

	    infile >> inittype;
	    infile >> initx;
	    if (geo == threeD) {
	      infile >> inity;
	      infile >> initz;
	    }
	    infile >> inite;
	    if (geo == threeD) {
	      infile >> inittheta;
	      infile >> initphi;
	    }
	    infile >> t;
	    if (geo == threeD) {
	      infile >> x;
	      infile >> y;
	      infile >> z;
	    }
	    infile >> energy;
	    if (geo == threeD) {
	      infile >> theta;
	      infile >> phi;
	    }

	    TLorentzVector v1(inite,0,0,inite);
	    v1.SetTheta(inittheta);
	    v1.SetPhi(initphi);
	    
	    TLorentzVector v2(energy,0,0,energy);
	    v2.SetTheta(theta);
	    v2.SetPhi(phi);
	    
	    ev.push_back( new UHECREvent((int)type, (int)inittype, TVector3(initx,inity,initz), v1, t, TVector3(x,y,z), v2) );
	  }
	  else infile.ignore(256,'\n');
	}
      }
    } 
    else Warning("CRPropaEvent", "Your ASCII file \"%s\" did not contain UHECR.", filename.c_str());

    infile.close();
    delete [] s;
  }
  else if (filetype == "FITS") {

    fitsfile* fitsinfile;
    int status = 0;
    string fitsfilename = filename;
    fitsfilename += "[EVENT_TABLE]";
    if (fits_open_file(&fitsinfile, filename.c_str(), 0, &status)) fits_report_error(stderr, status);
    if (fits_movnam_hdu(fitsinfile, BINARY_TBL, "EVENT_TABLE", 0, &status)) fits_report_error(stderr, status);
    else Ok = true;

    int ncols = 0;
    if (fits_get_num_cols(fitsinfile, &ncols, &status)) fits_report_error(stderr, status);
    if (ncols == 6) geo = oneD;            // CRPropa nuclei
    else if (ncols == 15) geo = threeD;    // CRPropa nuclei
    else Error("CRPropaEvent", "Wrong number of columns: %i", ncols);

    LONGLONG nrows = 0;
    if (fits_get_num_rowsll(fitsinfile, &nrows, &status)) fits_report_error(stderr, status);

    int* fitstype = new int[nrows]();
    int* fitsinittype = new int[nrows]();
    double* fitsinitposX = new double[nrows]();
    double* fitsinitposY = new double[nrows]();
    double* fitsinitposZ = new double[nrows]();
    double* fitsiniten = new double[nrows]();
    double* fitsinitentheta = new double[nrows]();
    double* fitsinitenphi = new double[nrows]();
    double* fitstime = new double[nrows]();
    double* fitsen = new double[nrows]();
    double* fitsentheta = new double[nrows]();
    double* fitsenphi = new double[nrows]();
    double* fitsposX = new double[nrows]();
    double* fitsposY = new double[nrows]();
    double* fitsposZ = new double[nrows]();
    LONGLONG firstrow = 1;
    LONGLONG firstelem = 1;
    int nulval = 0;
    double dnulval = 0;
    int anynul = 0;

    if (fits_read_col(fitsinfile, TINT, 1, firstrow, firstelem, nrows, &nulval, fitstype, &anynul, &status)) fits_report_error(stderr, status);
    if (fits_read_col(fitsinfile, TINT, 2, firstrow, firstelem, nrows, &nulval, fitsinittype, &anynul, &status)) fits_report_error(stderr, status);

    if (fits_read_col(fitsinfile, TDOUBLE, 3, firstrow, firstelem, nrows, &dnulval, fitsinitposX, &anynul, &status)) fits_report_error(stderr, status);

    if (geo == oneD) {
      if (fits_read_col(fitsinfile, TDOUBLE, 4, firstrow, firstelem, nrows, &dnulval, fitsiniten, &anynul, &status)) fits_report_error(stderr, status);

      if (fits_read_col(fitsinfile, TDOUBLE, 5, firstrow, firstelem, nrows, &dnulval, fitstime, &anynul, &status)) fits_report_error(stderr, status);

      if (fits_read_col(fitsinfile, TDOUBLE, 6, firstrow, firstelem, nrows, &dnulval, fitsen, &anynul, &status)) fits_report_error(stderr, status);
    }
    else {
      if (fits_read_col(fitsinfile, TDOUBLE, 4, firstrow, firstelem, nrows, &dnulval, fitsinitposY, &anynul, &status)) fits_report_error(stderr, status);
      if (fits_read_col(fitsinfile, TDOUBLE, 5, firstrow, firstelem, nrows, &dnulval, fitsinitposZ, &anynul, &status)) fits_report_error(stderr, status);

      if (fits_read_col(fitsinfile, TDOUBLE, 6, firstrow, firstelem, nrows, &dnulval, fitsiniten, &anynul, &status)) fits_report_error(stderr, status);
      if (fits_read_col(fitsinfile, TDOUBLE, 7, firstrow, firstelem, nrows, &dnulval, fitsinitentheta, &anynul, &status)) fits_report_error(stderr, status);
      if (fits_read_col(fitsinfile, TDOUBLE, 8, firstrow, firstelem, nrows, &dnulval, fitsinitenphi, &anynul, &status)) fits_report_error(stderr, status);

      if (fits_read_col(fitsinfile, TDOUBLE, 9, firstrow, firstelem, nrows, &dnulval, fitstime, &anynul, &status)) fits_report_error(stderr, status);

   if (fits_read_col(fitsinfile, TDOUBLE, 10, firstrow, firstelem, nrows, &dnulval, fitsposX, &anynul, &status)) fits_report_error(stderr, status);
   if (fits_read_col(fitsinfile, TDOUBLE, 11, firstrow, firstelem, nrows, &dnulval, fitsposY, &anynul, &status)) fits_report_error(stderr, status);
   if (fits_read_col(fitsinfile, TDOUBLE, 12, firstrow, firstelem, nrows, &dnulval, fitsposZ, &anynul, &status)) fits_report_error(stderr, status);

      if (fits_read_col(fitsinfile, TDOUBLE, 13, firstrow, firstelem, nrows, &dnulval, fitsen, &anynul, &status)) fits_report_error(stderr, status);
      if (fits_read_col(fitsinfile, TDOUBLE, 14, firstrow, firstelem, nrows, &dnulval, fitsentheta, &anynul, &status)) fits_report_error(stderr, status);
      if (fits_read_col(fitsinfile, TDOUBLE, 15, firstrow, firstelem, nrows, &dnulval, fitsenphi, &anynul, &status)) fits_report_error(stderr, status);
    }

    for (LONGLONG i = 0; i < nrows; ++i) {
      
      TLorentzVector v1(fitsiniten[i],0,0,fitsiniten[i]);
      v1.SetTheta(fitsinitentheta[i]);
      v1.SetPhi(fitsinitenphi[i]);
      
      TLorentzVector v2(fitsen[i],0,0,fitsen[i]);
      v2.SetTheta(fitsentheta[i]);
      v2.SetPhi(fitsenphi[i]);

      ev.push_back( new UHECREvent(fitstype[i], fitsinittype[i], TVector3(fitsinitposX[i],fitsinitposY[i],fitsinitposZ[i]), v1, fitstime[i], TVector3(fitsposX[i],fitsposY[i],fitsposZ[i]), v2) );
    }

    delete [] fitstype;
    delete [] fitsinittype;
    delete [] fitsinitposX;
    delete [] fitsinitposY;
    delete [] fitsinitposZ;
    delete [] fitsiniten;
    delete [] fitsinitentheta;
    delete [] fitsinitenphi;
    delete [] fitstime;
    delete [] fitsen;
    delete [] fitsentheta;
    delete [] fitsenphi;
    delete [] fitsposX;
    delete [] fitsposY;
    delete [] fitsposZ;

    if (fits_close_file(fitsinfile, &status)) fits_report_error(stderr, status);
  }
  else Fatal("CRPropaEvent",  "Unknown file type \"%s\" in CRPropaEvent.", filetype.c_str());

  whoami = "CR";

  return ;
}

TH1D* CRPropaEvent::ComputeSpectrum(const char* name, const double& newindex, const double& oldindex, const int& nbins, const double& Emin, const double& Emax, const double& factor, const double& reweightfactor) {
  
  if (Emax < Emin) {
    Error("ComputeSpectrum", "Wrong energy limits. Emin = %f. Emax = %f", Emin, Emax);
    return NULL;
  }

  TH1D* h1 = new TH1D(name, "", nbins, Emin, Emax);
  h1->Sumw2();
  for (unsigned int i = 0; i < ev.size(); ++i) {
    h1->Fill(18.0+log10(ev[i]->GetEnergy()),pow(1e18*ev[i]->GetEnergy(),factor)*pow(1.e18*ev[i]->GetInEnergy(),(oldindex-newindex)));
  }
  return h1;
}

TProfile* CRPropaEvent::ComputeMassSpectrum(const char* name, const double& newindex, const double& oldindex, const int& nbins, const double& Emin, const double& Emax) { 

  if (Emax < Emin) {
    Error("ComputeMassSpectrum", "Wrong energy limits. Emin = %f. Emax = %f", Emin, Emax);
    return NULL;
  }

  TProfile* h1 = new TProfile(name, "", nbins, Emin, Emax);
  h1->Sumw2();
  for (unsigned int i = 0; i < ev.size(); ++i) h1->Fill(log10(1e18*ev[i]->GetEnergy()), double(ev[i]->GetType()%100), pow(1.e18*ev[i]->GetInEnergy(),(oldindex-newindex)));
  return h1; 
}

TProfile* CRPropaEvent::ComputeDeflections(const char* name, const double& newindex, const double& oldindex, const int& nbins, const double& Emin, const double& Emax) {

  if (Emax < Emin) {
    Error("ComputeDeflections", "Wrong energy limits. Emin = %f. Emax = %f", Emin, Emax);
    return NULL;
  }

  TProfile* h1 = new TProfile(name, "", nbins, Emin, Emax);
  h1->Sumw2();
  for (unsigned int i = 0; i < ev.size(); ++i) {
    h1->Fill(18.0+log10(ev[i]->GetEnergy()), 180.0-ev[i]->GetMom().Angle(ev[i]->GetInPos())*TMath::RadToDeg(), pow(1.e18*ev[i]->GetInEnergy(),(oldindex-newindex)));
  }

  return h1;
}

TProfile* CRPropaEvent::ComputeCumulativeDeflections(const char* name, const double& newindex, const double& oldindex, const int& nbins, const double& Emin, const double& Emax, const double& Dmax) { 
  if (Emax < Emin) {
    Error("ComputeCumulativeDeflections", "Wrong energy limits. Emin = %f. Emax = %f", Emin, Emax);
    return NULL;
  }

  TProfile* h1 = new TProfile(name, "", nbins, 0, 180);
  h1->Sumw2();

  double angle = 0.0;

  for (int j = 0; j < nbins; ++j) {
    angle = double(j)*180.0/double(nbins-1);
    for (unsigned int i = 0; i < ev.size(); ++i) h1->Fill(angle, (ev[i]->GetEnergy() > Emin && ev[i]->GetEnergy() < Emax && ev[i]->GetTime() < Dmax && 180.0-ev[i]->GetMom().Angle(ev[i]->GetInPos())*TMath::RadToDeg() >= angle), pow(ev[i]->GetInEnergy(),(oldindex-newindex)));
  }

  return h1;
}

vector<double> CRPropaEvent::ComputeSkyMap(const double& newindex, const double& oldindex, const double& Emin, const double& Emax, const double& Dmax, const long& nside) {

  const long npix = nside2npix(nside);

  vector<double> map(npix, 0.0);

  for (unsigned int i = 0; i < ev.size(); ++i) {
    if (ev[i]->GetEnergy() > Emin && ev[i]->GetEnergy() < Emax && ev[i]->GetTime() < Dmax) {

      TVector3 v(-ev[i]->GetMom().Vect().Unit());

      double coord[3];
      coord[0] = v(0);
      coord[1] = v(1);
      coord[2] = v(2);

      long ipix = -1;
      vec2pix_ring(nside, coord, &ipix);

      if (ipix < 0 || ipix >= npix) Fatal("ComputeSkyMap", "ipix was out of bounds: ipix = %li", ipix);
      map[ipix] += pow(ev[i]->GetInEnergy(),(oldindex-newindex));
    }
  }
  return map;
}

void CRPropaEvent::SkyMapEvents(const char* name, const double& Emin, const double& Emax, const double& Dmax) {
  
  if (geo == oneD) {
    Fatal("SkyMapEvents", "We are in 1-D geometry. Your request does not make sense.");
    return ;
  }

  ofstream outfile(name, ios::out);
  
  for (unsigned int i = 0; i < ev.size(); ++i) {
    if (ev[i]->GetEnergy() > Emin && ev[i]->GetEnergy() < Emax && ev[i]->GetTime() < Dmax) {
      
      TVector3 v(-ev[i]->GetMom().Vect().Unit());

      double theta = v.Theta()*TMath::RadToDeg();
      double phi = v.Phi()* TMath::RadToDeg();

      double thetas = ev[i]->GetInPos().Theta()*TMath::RadToDeg();
      double phis = ev[i]->GetInPos().Phi()*TMath::RadToDeg();

      theta = 90-theta;
      phi = 180-phi;

      thetas = 90-thetas;
      phis = 180-phis;

      outfile <<  phi << " " << theta << " " << phis << " " << thetas << endl;
    }
  }
  outfile.close();
  return ;
}

void CRPropaEvent::SkyMapSources(const char* name, const double& Dmax) {

  if (geo == oneD) {
    Fatal("SkyMapSources", "We are in 1-D geometry. Your request does not make sense.");
    return ;
  }

  ofstream outfile(name, ios::out);

  vector<TVector3> InitPosition;
  for (unsigned int i = 0; i < ev.size(); ++i) InitPosition.push_back(ev[i]->GetInPos());

  for (vector<TVector3>::iterator i = InitPosition.begin(); i != InitPosition.end(); ++i) {
    if ((*i).Mag() < Dmax) {
      if ((int)count(InitPosition.begin(), i, (*i)) == 0) {

	double theta = (*i).Theta()*TMath::RadToDeg();
	double phi = (*i).Phi()* TMath::RadToDeg();
	if (phi < 0) phi += 360;
	outfile << phi << " " << theta <<endl;
      }
    }
  }
  outfile.close();
  return ;
}

/*****************************************************************
 *****************************************************************
 */

CRPropaNeutrino::CRPropaNeutrino(string& filename, const string& filetype) {

  whoami = "Neutrino";

  Ok = false;

  float type = 0;
  float inittype = 0;
  float initx = 0;
  float inity = 0;
  float initz = 0;
  float x = 0;
  float y = 0;
  float z = 0;
  float inite = 0;
  float energy = 0;
  float theta = 0;
  float phi = 0;
  float sourcenergy = 0;
  float sourceposx = 0;
  float sourceposy = 0;
  float sourceposz = 0;

  if (filename == "." || filename == "..") return ;
  if (filetype == "ROOT") {

    TFile rootfile(filename.c_str(), "READ");
    if (rootfile.IsOpen()) {
      TNtuple* t1;
      rootfile.GetObject("neutrinos",t1);

      if (t1) {

	Ok = true;
	TString title(t1->GetTitle());
	if (title.Contains("1D")) geo = oneD;
	else geo = threeD;

	t1->SetBranchAddress("Particle_Type", &type);
	t1->SetBranchAddress("Initial_Type", &inittype);
	t1->SetBranchAddress("Source_Energy_EeV", &sourcenergy);

	if (geo == oneD) {
	  t1->SetBranchAddress("Initial_Position_Mpc", &initx);
	  t1->SetBranchAddress("Energy_EeV", &energy);
	  t1->SetBranchAddress("Initial_Energy_EeV", &inite);
	  t1->SetBranchAddress("Source_Position_Mpc", &sourceposx);
	}
	else {
	  t1->SetBranchAddress("Initial_Position_X_Mpc", &initx);
	  t1->SetBranchAddress("Initial_Position_Y_Mpc", &inity);
	  t1->SetBranchAddress("Initial_Position_Z_Mpc", &initz);

	  t1->SetBranchAddress("Position_X_Mpc", &x);
	  t1->SetBranchAddress("Position_Y_Mpc", &y);
	  t1->SetBranchAddress("Position_Z_Mpc", &z);

	  t1->SetBranchAddress("Source_Position_X_Mpc", &sourceposx);
	  t1->SetBranchAddress("Source_Position_Y_Mpc", &sourceposy);
	  t1->SetBranchAddress("Source_Position_Z_Mpc", &sourceposz);

	  t1->SetBranchAddress("Momentum_E_EeV", &energy);
	  t1->SetBranchAddress("Momentum_theta", &theta);
	  t1->SetBranchAddress("Momentum_phi", &phi);
	}

	int Nev = t1->GetEntries();

	for (int i = 0; i < Nev; i++) {
	  t1->GetEntry(i);
	  TLorentzVector v2(energy,0,0,energy);
	  v2.SetTheta(theta);
	  v2.SetPhi(phi);

	  neu.push_back( new UHECRNeutrino((int)type, (int)inittype, TVector3(sourceposx,sourceposy,sourceposz), TVector3(initx,inity,initz), sourcenergy, TVector3(x,y,z), v2, sourcenergy) );
	}
	rootfile.Close();
	t1 = NULL;
      } // t1
      else Warning("CRPropaNeutrino", "Your ROOT file\"%s\" did not contain any neutrino.", filename.c_str());
    } // rootfile
    else Error("CRPropaNeutrino", "Could not open the ROOT file \"%s\".", filename.c_str());
  } // strcmp
  else if (filetype == "ASCII") {

    ifstream infile(filename.c_str(), ios::in);
    if (!infile.is_open()) Error("CRPropaNeutrino", "Could not open the ASCII file \"%s\".", filename.c_str());
    
    string str("");
    char* s = new char[3000]();
    
    while (strcmp(s, "#Output secondary neutrinos") && infile.good()) infile.getline(s,3000);
    if (infile.good()) {
      infile.getline(s,3000);
      TString rootstring(s);
      if (rootstring.Contains("[")) geo = threeD;
      else geo = oneD;
    
      while (infile.good()) {
	infile >> str;
	if (str[0] == '#') infile.ignore(256,'\n');
	else {
	  type = atoi(str.c_str());
	  if (abs(type) < 20) {
	    Ok = true;
	    infile >> inittype;
	    infile >> sourceposx;
	    if (geo == threeD) {
	      infile >> sourceposy;
	      infile >> sourceposz;
	    }
	    infile >> initx;
	    if (geo == threeD) {
	      infile >> inity;
	      infile >> initz;
	    }
	    if (geo == threeD) {
	      infile >> x;
	      infile >> y;
	      infile >> z;
	    }
	    else infile >> inite;
	    infile >> energy;
	    if (geo == threeD) {
	      infile >> theta;
	      infile >> phi;
	    }
	    infile >> sourcenergy;

	    TLorentzVector v2(energy,0,0,energy);
	    v2.SetTheta(theta);
	    v2.SetPhi(phi);

	    neu.push_back( new UHECRNeutrino((int)type, (int)inittype, TVector3(sourceposx,sourceposy,sourceposz), TVector3(initx,inity,initz), inite, TVector3(x,y,z), v2, sourcenergy) );
	  }
	}
      }
    } 
    else Warning("CRPropaNeutrino", "Your ASCII file\"%s\" did not contain any neutrino.", filename.c_str());

    infile.close();
    delete [] s;
  }
  else if (filetype == "FITS") {

    fitsfile* fitsinfile;
    int status = 0;
    string fitsfilename = filename;
    fitsfilename += "[NEUTRINO_TABLE]";
    if (fits_open_file(&fitsinfile, fitsfilename.c_str(), 0, &status)) fits_report_error(stderr, status);
    else Ok = true;

    int ncols = 0;
    if (fits_get_num_cols(fitsinfile, &ncols, &status)) fits_report_error(stderr, status);
    if (ncols == 7) geo = oneD;
    else if (ncols == 15) geo = threeD;
    else Error("CRPropaNeutrino", "Wrong number of columns: %i", ncols);

    LONGLONG nrows = 0;
    if (fits_get_num_rowsll(fitsinfile, &nrows, &status)) fits_report_error(stderr, status);

    int* fitstype = new int[nrows]();
    int* fitsinittype = new int[nrows]();
    double* fitssourceposX = new double[nrows]();
    double* fitssourceposY = new double[nrows]();
    double* fitssourceposZ = new double[nrows]();
    double* fitsinitposX = new double[nrows]();
    double* fitsinitposY = new double[nrows]();
    double* fitsinitposZ = new double[nrows]();
    double* fitsiniten = new double[nrows]();
    double* fitssourceen = new double[nrows]();
    double* fitsen = new double[nrows]();
    double* fitsentheta = new double[nrows]();
    double* fitsenphi = new double[nrows]();
    double* fitsposX = new double[nrows]();
    double* fitsposY = new double[nrows]();
    double* fitsposZ = new double[nrows]();
    LONGLONG firstrow = 1;
    LONGLONG firstelem = 1;
    int nulval = 0;
    double dnulval = 0;
    int anynul = 0;
    
    if (fits_read_col(fitsinfile, TINT, 1, firstrow, firstelem, nrows, &nulval, fitstype, &anynul, &status)) fits_report_error(stderr, status);
    if (fits_read_col(fitsinfile, TINT, 2, firstrow, firstelem, nrows, &nulval, fitsinittype, &anynul, &status)) fits_report_error(stderr, status);

    if (fits_read_col(fitsinfile, TDOUBLE, 3, firstrow, firstelem, nrows, &dnulval, fitssourceposX, &anynul, &status)) fits_report_error(stderr, status);

    if (geo == oneD) {
      if (fits_read_col(fitsinfile, TDOUBLE, 4, firstrow, firstelem, nrows, &dnulval, fitsinitposX, &anynul, &status)) fits_report_error(stderr, status);

      if (fits_read_col(fitsinfile, TDOUBLE, 5, firstrow, firstelem, nrows, &dnulval, fitsiniten, &anynul, &status)) fits_report_error(stderr, status);

      if (fits_read_col(fitsinfile, TDOUBLE, 6, firstrow, firstelem, nrows, &dnulval, fitsen, &anynul, &status)) fits_report_error(stderr, status);

      if (fits_read_col(fitsinfile, TDOUBLE, 7, firstrow, firstelem, nrows, &dnulval, fitssourceen, &anynul, &status)) fits_report_error(stderr, status);
    }
    else {
      if (fits_read_col(fitsinfile, TDOUBLE, 4, firstrow, firstelem, nrows, &dnulval, fitssourceposY, &anynul, &status)) fits_report_error(stderr, status);
      if (fits_read_col(fitsinfile, TDOUBLE, 5, firstrow, firstelem, nrows, &dnulval, fitssourceposZ, &anynul, &status)) fits_report_error(stderr, status);

      if (fits_read_col(fitsinfile, TDOUBLE, 6, firstrow, firstelem, nrows, &dnulval, fitsinitposX, &anynul, &status)) fits_report_error(stderr, status);
      if (fits_read_col(fitsinfile, TDOUBLE, 7, firstrow, firstelem, nrows, &dnulval, fitsinitposY, &anynul, &status)) fits_report_error(stderr, status);
      if (fits_read_col(fitsinfile, TDOUBLE, 8, firstrow, firstelem, nrows, &dnulval, fitsinitposZ, &anynul, &status)) fits_report_error(stderr, status);

   if (fits_read_col(fitsinfile, TDOUBLE, 9, firstrow, firstelem, nrows, &dnulval, fitsposX, &anynul, &status)) fits_report_error(stderr, status);
   if (fits_read_col(fitsinfile, TDOUBLE, 10, firstrow, firstelem, nrows, &dnulval, fitsposY, &anynul, &status)) fits_report_error(stderr, status);
   if (fits_read_col(fitsinfile, TDOUBLE, 11, firstrow, firstelem, nrows, &dnulval, fitsposZ, &anynul, &status)) fits_report_error(stderr, status);

      if (fits_read_col(fitsinfile, TDOUBLE, 12, firstrow, firstelem, nrows, &dnulval, fitsen, &anynul, &status)) fits_report_error(stderr, status);
      if (fits_read_col(fitsinfile, TDOUBLE, 13, firstrow, firstelem, nrows, &dnulval, fitsentheta, &anynul, &status)) fits_report_error(stderr, status);
      if (fits_read_col(fitsinfile, TDOUBLE, 14, firstrow, firstelem, nrows, &dnulval, fitsenphi, &anynul, &status)) fits_report_error(stderr, status);

      if (fits_read_col(fitsinfile, TDOUBLE, 14, firstrow, firstelem, nrows, &dnulval, fitssourceen, &anynul, &status)) fits_report_error(stderr, status);
    }

    for (LONGLONG i = 0; i < nrows; ++i) {
      
      TLorentzVector v2(fitsen[i],0,0,fitsen[i]);
      v2.SetTheta(fitsentheta[i]);
      v2.SetPhi(fitsenphi[i]);

      neu.push_back( new UHECRNeutrino(fitstype[i], fitsinittype[i], TVector3(fitssourceposX[i],fitssourceposY[i],fitssourceposZ[i]), TVector3(fitsinitposX[i],fitsinitposY[i],fitsinitposZ[i]), fitsiniten[i], TVector3(fitsposX[i],fitsposY[i],fitsposZ[i]), v2, fitssourceen[i]) );
    }

    delete [] fitstype;
    delete [] fitsinittype;
    delete [] fitssourceposX;
    delete [] fitssourceposY;
    delete [] fitssourceposZ;
    delete [] fitsinitposX;
    delete [] fitsinitposY;
    delete [] fitsinitposZ;
    delete [] fitsiniten;
    delete [] fitssourceen;
    delete [] fitsen;
    delete [] fitsentheta;
    delete [] fitsenphi;
    delete [] fitsposX;
    delete [] fitsposY;
    delete [] fitsposZ;

    if (fits_close_file(fitsinfile, &status)) fits_report_error(stderr, status);
  }
  else Fatal("CRPropaNeutrino",  "Unknown file type \"%s\".", filetype.c_str());

  return ;
}

TH1D* CRPropaNeutrino::ComputeSpectrum(const char* name, const double& newindex, const double& oldindex, const int& nbins, const double& Emin, const double& Emax, const double& factor, const double& reweightfactor) {
  if (Emax < Emin) {
    Error("ComputeSpectrum", "Wrong energy limits. Emin = %f. Emax = %f", Emin, Emax);
    return NULL;
  }
 
  TH1D* h1 = new TH1D(name, "", nbins, Emin, Emax);
  h1->Sumw2();
  for (unsigned int i = 0; i < neu.size(); ++i) h1->Fill(18.0+log10(neu[i]->GetEnergy()), pow(1.e18*neu[i]->GetSourceEnergy(),(oldindex-newindex))*pow(1.e18*neu[i]->GetEnergy(), factor));

  return h1;
}

/*****************************************************************
 *****************************************************************
 */

CRPropaPhotonSpectrum::CRPropaPhotonSpectrum(string& filename, const string& filetype) {

  whoami = "PhotonSpectrum";

  Ok = false;

  int type = 0;
  int inittype = 0;
  double initx = 0;
  double inity = 0;
  double initz = 0;
  double x = 0;
  double y = 0;
  double z = 0;
  char* cascorig = NULL;
  double energy = 0;
  double theta = 0;
  double phi = 0;
  double sourcenergy = 0;
  double sourceposx = 0;
  double sourceposy = 0;
  double sourceposz = 0;
  float spectrum[170];

  if (filename == "." || filename == "..") return ;
  if (filetype == "ROOT") {

    TFile rootfile(filename.c_str(), "READ");
    if (rootfile.IsOpen()) {
      TTree* t1;
      rootfile.GetObject("photons",t1);

      if (t1) {

	char* cascorigROOT = new char[20];//NULL;
	Ok = true;

	TString title(t1->GetTitle());
	if (title.Contains("1D")) geo = oneD;
	else geo = threeD;

	t1->SetBranchAddress("Particle_Type", &type);
	t1->SetBranchAddress("Initial_Type", &inittype);
	t1->SetBranchAddress("Cascade_Origin_String", cascorigROOT);
	t1->SetBranchAddress("SourceEnergy_EeV", &sourcenergy);
	t1->SetBranchAddress("Spectrum", &spectrum);
	  
	if (geo == oneD) {
	  t1->SetBranchAddress("InitialPosition_Mpc", &initx);
	  t1->SetBranchAddress("SourcePosition_Mpc", &sourceposx);
	}
	else {
	  t1->SetBranchAddress("InitialPosition_X_Mpc", &initx);
	  t1->SetBranchAddress("InitialPosition_Y_Mpc", &inity);
	  t1->SetBranchAddress("InitialPosition_Z_Mpc", &initz);

	  t1->SetBranchAddress("Position_X_Mpc", &x);
	  t1->SetBranchAddress("Position_Y_Mpc", &y);
	  t1->SetBranchAddress("Position_Z_Mpc", &z);

	  t1->SetBranchAddress("SourcePosition_X_Mpc", &sourceposx);
	  t1->SetBranchAddress("SourcePosition_Y_Mpc", &sourceposy);
	  t1->SetBranchAddress("SourcePosition_Z_Mpc", &sourceposz);

	  t1->SetBranchAddress("Energy_EeV", &energy);
	  t1->SetBranchAddress("Momentum_theta", &theta);
	  t1->SetBranchAddress("Momentum_phi", &phi);
	}

	int Nev = t1->GetEntries();

	for (int i = 0; i < Nev; i++) {
	  t1->GetEntry(i);

	  TLorentzVector v2(energy,0,0,energy);
	  v2.SetTheta(theta);
	  v2.SetPhi(phi);

	  vector<double> result;
	  for (int j = 0; j < 170; j++) result.push_back(spectrum[j]);

	  phot.push_back( new UHECRPhoton(type, inittype, std::string(cascorigROOT), TVector3(sourceposx,sourceposy,sourceposz), TVector3(initx,inity,initz), TVector3(x,y,z), v2, sourcenergy, result) );
	}
	rootfile.Close();
	t1 = NULL;
	delete [] cascorigROOT;
	cascorigROOT=NULL;
      } // t1
      else Warning("CRPropaPhotonSpectrum",  "Your ROOT file \"%s\" did not contain photon showers.", filename.c_str());
    } // rootfile
    else Error("CRPropaPhotonSpectrum",  "Could not open your ROOT file \"%s\".", filename.c_str());
  } // strcmp
  else if (filetype == "ASCII") {

    ifstream infile(filename.c_str(), ios::in);
    if (!infile.is_open()) Error("CRPropaPhotonSpectrum",  "Could not open your ASCII file \"%s\".", filename.c_str());

    string str("");
    char* s = new char[3000]();
    
    while (strcmp(s, "#Output secondary photons") && infile.good()) infile.getline(s,3000);
    if (infile.good()) {
      infile.getline(s,3000);
      TString rootstring(s);
      if (rootstring.Contains("[")) geo = threeD;
      else geo = oneD;

      while (infile.good()) {
	infile >> str;
	if (str[0] == '#') {
	  infile >> str; infile >> str;
	  if (str.find("neutrino") != string::npos) break;
	  infile.ignore(516,'\n');
	}
	else {
	  Ok = true;

	  type = atoi(str.c_str());
	  if (type == 22) {
	    infile >> str;
	    if (!strcmp(str.c_str(), "Ebins_EeV")) {
	      infile.getline(s,3000);
	      infile.getline(s,3000);
	    }
	    else {

	      cascorig = (char*)str.c_str();
	      infile >> inittype;

	      infile >> sourceposx;
	      if (geo == threeD) {
		infile >> sourceposy;
		infile >> sourceposz;
	      }
	      infile >> initx;
	      if (geo == threeD) {
		infile >> inity;
		infile >> initz;
	      }
	      if (geo == threeD) {
		infile >> x;
		infile >> y;
		infile >> z;
	      }
	      infile >> energy;
	      if (geo == threeD) {
		infile >> theta;
		infile >> phi;
	      }
	      infile >> sourcenergy;

	      vector<double> result(170,0.0);
	      for (int i = 0; i < 170; i++) infile >> result[i];

	      TLorentzVector v2(energy,0,0,energy);
	      v2.SetTheta(theta);
	      v2.SetPhi(phi);
	    
	      phot.push_back( new UHECRPhoton(type, inittype, std::string(cascorig), TVector3(sourceposx,sourceposy,sourceposz), TVector3(initx,inity,initz), TVector3(x,y,z), v2, sourcenergy, result) );

	    }
	  }
	}
      }
    } 
    else Warning("CRPropaPhotonSpectrum",  "Your ASCII file \"%s\" did not contain photon showers.", filename.c_str());

    infile.close();
    delete [] s;
  }
  else if (filetype == "FITS") {

    fitsfile* fitsinfile;
    int status = 0;
    string fitsfilename = filename;
    fitsfilename += "[PHOTON_SPECTRUM_TABLE]";
    if (fits_open_file(&fitsinfile, fitsfilename.c_str(), 0, &status)) fits_report_error(stderr, status);
    else Ok = true;

    int ncols = 0;
    if (fits_get_num_cols(fitsinfile, &ncols, &status)) fits_report_error(stderr, status);
    if (ncols == 7) geo = oneD;
    else if (ncols == 17) geo = threeD;
    else Error("CRPropaPhotonSpectrum", "Wrong number of columns: %i", ncols);

    LONGLONG nrows = 0;
    if (fits_get_num_rowsll(fitsinfile, &nrows, &status)) fits_report_error(stderr, status);

    int* fitstype = new int[nrows]();
    int* fitsinittype = new int[nrows]();
    char** fitscascorig = new char*[nrows]();
    for (LONGLONG i = 0; i < nrows; ++i) fitscascorig[i] = new char[11]();
    double* fitssourceposX = new double[nrows]();
    double* fitssourceposY = new double[nrows]();
    double* fitssourceposZ = new double[nrows]();
    double* fitsinitposX = new double[nrows]();
    double* fitsinitposY = new double[nrows]();
    double* fitsinitposZ = new double[nrows]();
    double* fitssourceen = new double[nrows]();
    double* fitsen = new double[nrows]();
    double* fitsentheta = new double[nrows]();
    double* fitsenphi = new double[nrows]();
    double* fitsposX = new double[nrows]();
    double* fitsposY = new double[nrows]();
    double* fitsposZ = new double[nrows]();
    double* fitsspectrum = new double[170*nrows]();

    LONGLONG firstrow = 1;
    LONGLONG firstelem = 1;
    int nulval = 0;
    double dnulval = 0;
    int anynul = 0;

    if (fits_read_col(fitsinfile, TINT, 1, firstrow, firstelem, nrows, &nulval, fitstype, &anynul, &status)) fits_report_error(stderr, status);
    if (fits_read_col(fitsinfile, TSTRING, 2, firstrow, firstelem, nrows, &nulval, fitscascorig, &anynul, &status)) fits_report_error(stderr, status);

    if (geo == oneD) {
      if (fits_read_col(fitsinfile, TDOUBLE, 3, firstrow, firstelem, nrows, &dnulval, fitssourceposX, &anynul, &status)) fits_report_error(stderr, status);
      if (fits_read_col(fitsinfile, TDOUBLE, 4, firstrow, firstelem, nrows, &dnulval, fitsinitposX, &anynul, &status)) fits_report_error(stderr, status);

      if (fits_read_col(fitsinfile, TDOUBLE, 5, firstrow, firstelem, nrows, &dnulval, fitsen, &anynul, &status)) fits_report_error(stderr, status);

      if (fits_read_col(fitsinfile, TDOUBLE, 6, firstrow, firstelem, nrows, &dnulval, fitssourceen, &anynul, &status)) fits_report_error(stderr, status);
    }
    else {
      if (fits_read_col(fitsinfile, TINT, 3, firstrow, firstelem, nrows, &nulval, fitsinittype, &anynul, &status)) fits_report_error(stderr, status);

      if (fits_read_col(fitsinfile, TDOUBLE, 4, firstrow, firstelem, nrows, &dnulval, fitssourceposX, &anynul, &status)) fits_report_error(stderr, status);
      if (fits_read_col(fitsinfile, TDOUBLE, 5, firstrow, firstelem, nrows, &dnulval, fitssourceposY, &anynul, &status)) fits_report_error(stderr, status);
      if (fits_read_col(fitsinfile, TDOUBLE, 6, firstrow, firstelem, nrows, &dnulval, fitssourceposZ, &anynul, &status)) fits_report_error(stderr, status);

      if (fits_read_col(fitsinfile, TDOUBLE, 7, firstrow, firstelem, nrows, &dnulval, fitsinitposX, &anynul, &status)) fits_report_error(stderr, status);
      if (fits_read_col(fitsinfile, TDOUBLE, 8, firstrow, firstelem, nrows, &dnulval, fitsinitposY, &anynul, &status)) fits_report_error(stderr, status);
      if (fits_read_col(fitsinfile, TDOUBLE, 9, firstrow, firstelem, nrows, &dnulval, fitsinitposZ, &anynul, &status)) fits_report_error(stderr, status);

      if (fits_read_col(fitsinfile, TDOUBLE, 10, firstrow, firstelem, nrows, &dnulval, fitsposX, &anynul, &status)) fits_report_error(stderr, status);
      if (fits_read_col(fitsinfile, TDOUBLE, 11, firstrow, firstelem, nrows, &dnulval, fitsposY, &anynul, &status)) fits_report_error(stderr, status);
      if (fits_read_col(fitsinfile, TDOUBLE, 12, firstrow, firstelem, nrows, &dnulval, fitsposZ, &anynul, &status)) fits_report_error(stderr, status);

      if (fits_read_col(fitsinfile, TDOUBLE, 13, firstrow, firstelem, nrows, &dnulval, fitsen, &anynul, &status)) fits_report_error(stderr, status);
      if (fits_read_col(fitsinfile, TDOUBLE, 14, firstrow, firstelem, nrows, &dnulval, fitsentheta, &anynul, &status)) fits_report_error(stderr, status);
      if (fits_read_col(fitsinfile, TDOUBLE, 15, firstrow, firstelem, nrows, &dnulval, fitsenphi, &anynul, &status)) fits_report_error(stderr, status);

      if (fits_read_col(fitsinfile, TDOUBLE, 16, firstrow, firstelem, nrows, &dnulval, fitssourceen, &anynul, &status)) fits_report_error(stderr, status);
 
     }
      if (fits_read_col(fitsinfile, TDOUBLE, ncols, firstrow, firstelem, 170*nrows, &dnulval, fitsspectrum, &anynul, &status)) fits_report_error(stderr, status);

    for (LONGLONG i = 1; i < nrows; i++) {
      
      vector<double> result(170,0.0);
      for (LONGLONG j = 0; j < 170; j++) result[j] = fitsspectrum[i*170+j];

      TLorentzVector v2(fitsen[i],0,0,fitsen[i]);
      v2.SetTheta(fitsentheta[i]);
      v2.SetPhi(fitsenphi[i]);

      std::string strcascorig = fitscascorig[i];
      phot.push_back( new UHECRPhoton(fitstype[i], fitsinittype[i], strcascorig, TVector3(fitssourceposX[i],fitssourceposY[i],fitssourceposZ[i]), TVector3(fitsinitposX[i],fitsinitposY[i],fitsinitposZ[i]), TVector3(fitsposX[i],fitsposY[i],fitsposZ[i]), v2, fitssourceen[i], result) );
    }

    delete [] fitstype;
    for (LONGLONG i = 0; i < nrows; ++i) delete [] fitscascorig[i];
    delete [] fitscascorig;
    delete [] fitsinittype;
    delete [] fitssourceposX;
    delete [] fitssourceposY;
    delete [] fitssourceposZ;
    delete [] fitsinitposX;
    delete [] fitsinitposY;
    delete [] fitsinitposZ;
    delete [] fitssourceen;
    delete [] fitsen;
    delete [] fitsentheta;
    delete [] fitsenphi;
    delete [] fitsposX;
    delete [] fitsposY;
    delete [] fitsposZ;
    delete [] fitsspectrum;

    if (fits_close_file(fitsinfile, &status)) fits_report_error(stderr, status);
  }
  else Fatal("CRPropaPhotonSpectrum",  "Unknown file type \"%s\".", filetype.c_str());

  return ;
}

TH1D* CRPropaPhotonSpectrum::ComputeSpectrum(const char* name, const double& newindex, const double& oldindex, const int& nbins, const double& Emin, const double& Emax, const double& factor, const double& reweightfactor) {
  if (Emax < Emin) {
    Error("ComputeSpectrum", "Wrong energy limits. Emin = %f. Emax = %f", Emin, Emax);
    return NULL;
  }

  if (geo == oneD && fabs(newindex-oldindex) >= 0.1) Warning("ComputeSpectrum", "Reweighting of photon spectra in 1D is not possible. The original spectrum will be produced.");

  TH1D* h1 = new TH1D(name, "", nbins, Emin, Emax);
  h1->Sumw2();
  for (unsigned int i = 0; i < phot.size(); ++i) {
    for (unsigned int j = 0; j < 170; ++j) {
      if (geo == oneD) h1->Fill(18.0+log10(EnergyBins[j]), phot[i]->GetSpectrum()[j]*pow(1.e18*EnergyBins[j],factor));
      else {
	if (phot[i]->GetCascadeOrig() == "PairProd") h1->Fill(18.0+log10(EnergyBins[j]), reweightfactor*phot[i]->GetSpectrum()[j]*pow(1.e18*EnergyBins[j],factor)*pow(1.e18*phot[i]->GetSourceEnergy(),(oldindex-newindex)));
	else h1->Fill(18.0+log10(EnergyBins[j]), phot[i]->GetSpectrum()[j]*pow(1.e18*EnergyBins[j],factor)*pow(1.e18*phot[i]->GetSourceEnergy(),(oldindex-newindex)));
	//	cout << phot[i]->GetCascadeOrig() << " " << 18.0+log10(EnergyBins[j]) << " " << phot[i]->GetSpectrum()[j] << " " << pow(1.e18*EnergyBins[j],factor) << " " << 1.e18*phot[i]->GetSourceEnergy() << " " << pow(1.e18*phot[i]->GetSourceEnergy(),(oldindex-newindex)) << endl;
      }
    }
  }

  return h1;
}
/*
void CRPropaPhotonSpectrum::SkyMapEvents(const char* name, const double& Emin, const double& Emax, const double& Dmax) {

  if (geo == oneD) {
    Fatal("SkyMapEvents", "We are in 1-D geometry. Your request does not make sense.");
    return ;
  }

  ofstream outfile(name, ios::out);
  
  for (unsigned int i = 0; i < Momentum.size(); ++i) {
      
    TVector3 v(-Momentum[i].Vect().Unit());
    
    double theta = v.Theta()*TMath::RadToDeg();
    double phi = v.Phi()* TMath::RadToDeg();
    
    double thetas = SourcePosition[i].Theta()*TMath::RadToDeg();
    double phis = SourcePosition[i].Phi()*TMath::RadToDeg();
    
    theta = 90-theta;
    phi = 180-phi;
    
    thetas = 90-thetas;
    phis = 180-phis;
    
    outfile <<  phi << " " << theta << " " << phis << " " << thetas << endl;
  }
  outfile.close();
  return ;
}
*/
void CRPropaPhotonSpectrum::SkyMapSources(const char* name, const double& Dmax) {

  if (geo == oneD) {
    Fatal("SkyMapSources", "We are in 1-D geometry. Your request does not make sense.");
    return ;
  }

  ofstream outfile(name, ios::out);

  vector<TVector3> SourcePosition;
  for (unsigned int i = 0; i < phot.size(); ++i) SourcePosition.push_back(phot[i]->GetSourcePos());

  for (vector<TVector3>::iterator i = SourcePosition.begin(); i != SourcePosition.end(); ++i) {
    if ((*i).Mag() < Dmax) {
      if ((int)count(SourcePosition.begin(), i, (*i)) == 0) {

	double theta = (*i).Theta()*TMath::RadToDeg();
	double phi = (*i).Phi()* TMath::RadToDeg();
	if (phi < 0) phi += 360;
	outfile << phi << " " << theta <<endl;
      }
    }
  }
  outfile.close();
  return ;
}

/*****************************************************************
 *****************************************************************
 */

CRPropaTrajectory::CRPropaTrajectory(string& filename, const string& filetype) {

  whoami = "Trajectory";

  Ok = false;

  float type = 0;
  float inittype = 0;
  double x = 0;
  double y = 0;
  double z = 0;
  double energy = 0;
  double t = 0;

  if (filetype == "ROOT") {

    TFile rootfile(filename.c_str(), "READ");
    if (rootfile.IsOpen()) {
      TNtuple* t1;
      rootfile.GetObject("traj",t1);

      if (t1) {

	Ok = true;

	TString title(t1->GetTitle());
	if (title.Contains("1D")) geo = oneD;
	else geo = threeD;

	t1->SetBranchAddress("Particle_Type", &type);
	t1->SetBranchAddress("Initial_Type", &inittype);
	t1->SetBranchAddress("Time_Mpc", &t);
	t1->SetBranchAddress("Energy_EeV", &energy);

	  t1->SetBranchAddress("Position_X_Mpc", &x);
	if (geo == threeD) {
	  t1->SetBranchAddress("Position_Y_Mpc", &y);
	  t1->SetBranchAddress("Position_Z_Mpc", &z);
	}

	int Nev = t1->GetEntries();

	for (int i = 0; i < Nev; i++) {
	  t1->GetEntry(i);
	  if (type == -1 && inittype == -1) ntracks.push_back(i);
	  traj.push_back( new UHECRTrajectory((int)type, (int)inittype, t, TVector3(x,y,z), energy) );
	}

	rootfile.Close();
	t1 = NULL;
      } // t1
      else Warning("CRPropaTrajectory", "Your ROOT file\"%s\" did not contain any trajectory.", filename.c_str());
    } // rootfile
    else Error("CRPropaTrajectory", "Could not open your ROOT file\"%s\".", filename.c_str());
  } // strcmp
  else if (filetype == "ASCII") {

    ifstream infile(filename.c_str(), ios::in);
    if (!infile.is_open()) Error("CRPropaTrajectory", "Could not open your ASCII file\"%s\".", filename.c_str());
    
    string str("");
    char* s = new char[3000]();
    
    while (strcmp(str.c_str(), "#Format") && infile.good()) infile >> str;
    if (infile.good()) {
      infile.getline(s,3000);
      TString rootstring(s);
      if (rootstring.Contains("[")) geo = threeD;
      else geo = oneD;
    
      while (infile.good()) {
	infile >> str;
	if (str[0] == '#') infile.ignore(256,'\n');
	else {
	  type = atoi(str.c_str());
	  Ok = true;
	  infile >> inittype;
	  infile >> t;
	  infile >> x;
	  if (geo == threeD) {
	    infile >> y;
	    infile >> z;
	  }
	  infile >> energy;
	    
	  traj.push_back( new UHECRTrajectory((int)type, (int)inittype, t, TVector3(x,y,z), energy) );
	  if (type == -1 && inittype == -1) ntracks.push_back(traj.size()-1);
	}
      }
    }
    else Warning("CRPropaTrajectory", "Your ASCII file\"%s\" did not contain any trajectory.", filename.c_str());

    infile.close();
    delete [] s;
  }
  else if (filetype == "FITS") {

    fitsfile* fitsinfile;
    int status = 0;
    string fitsfilename = filename;
    fitsfilename += "[TRAJECTORY_TABLE]";
    if (fits_open_file(&fitsinfile, fitsfilename.c_str(), 0, &status)) fits_report_error(stderr, status);
    else Ok = true;

    int ncols = 0;
    if (fits_get_num_cols(fitsinfile, &ncols, &status)) fits_report_error(stderr, status);
    if (ncols == 5) geo = oneD;
    else if (ncols == 7) geo = threeD;
    else Error("CRPropaTrajectory", "Wrong number of columns: %i", ncols);

    LONGLONG nrows = 0;
    if (fits_get_num_rowsll(fitsinfile, &nrows, &status)) fits_report_error(stderr, status);

    int* fitstype = new int[nrows]();
    int* fitsinittype = new int[nrows]();
    double* fitstime = new double[nrows]();
    double* fitsen = new double[nrows]();
    double* fitsposX = new double[nrows]();
    double* fitsposY = new double[nrows]();
    double* fitsposZ = new double[nrows]();
    LONGLONG firstrow = 1;
    LONGLONG firstelem = 1;
    int nulval = 0;
    double dnulval = 0;
    int anynul = 0;

    if (fits_read_col(fitsinfile, TINT, 1, firstrow, firstelem, nrows, &nulval, fitstype, &anynul, &status)) fits_report_error(stderr, status);
    if (fits_read_col(fitsinfile, TINT, 2, firstrow, firstelem, nrows, &nulval, fitsinittype, &anynul, &status)) fits_report_error(stderr, status);

 
    if (geo == oneD) {
      if (fits_read_col(fitsinfile, TDOUBLE, 3, firstrow, firstelem, nrows, &dnulval, fitstime, &anynul, &status)) fits_report_error(stderr, status);
      if (fits_read_col(fitsinfile, TDOUBLE, 4, firstrow, firstelem, nrows, &dnulval, fitsposX, &anynul, &status)) fits_report_error(stderr, status);
      if (fits_read_col(fitsinfile, TDOUBLE, 5, firstrow, firstelem, nrows, &dnulval, fitsen, &anynul, &status)) fits_report_error(stderr, status);
    }
    else {
      if (fits_read_col(fitsinfile, TDOUBLE, 3, firstrow, firstelem, nrows, &dnulval, fitstime, &anynul, &status)) fits_report_error(stderr, status);

   if (fits_read_col(fitsinfile, TDOUBLE, 4, firstrow, firstelem, nrows, &dnulval, fitsposX, &anynul, &status)) fits_report_error(stderr, status);
   if (fits_read_col(fitsinfile, TDOUBLE, 5, firstrow, firstelem, nrows, &dnulval, fitsposY, &anynul, &status)) fits_report_error(stderr, status);
   if (fits_read_col(fitsinfile, TDOUBLE, 6, firstrow, firstelem, nrows, &dnulval, fitsposZ, &anynul, &status)) fits_report_error(stderr, status);

      if (fits_read_col(fitsinfile, TDOUBLE, 7, firstrow, firstelem, nrows, &dnulval, fitsen, &anynul, &status)) fits_report_error(stderr, status);
    }

    for (LONGLONG i = 0; i < nrows; ++i) {
      traj.push_back( new UHECRTrajectory(fitstype[i], fitsinittype[i], fitstime[i], TVector3(fitsposX[i],fitsposY[i],fitsposZ[i]), fitsen[i]) );
      if (type == -1 && inittype == -1) ntracks.push_back(traj.size()-1);
    }

    delete [] fitstype;
    delete [] fitsinittype;
    delete [] fitstime;
    delete [] fitsen;
    delete [] fitsposX;
    delete [] fitsposY;
    delete [] fitsposZ;

    if (fits_close_file(fitsinfile, &status)) fits_report_error(stderr, status);
  }
  else Fatal("CRPropaTrajectory", "Unknown file type \"%s\".", filetype.c_str());

  return ;
}

void CRPropaTrajectory::DrawTrajectory(unsigned int id) {

  
  if (id == 0) Fatal("DrawTrajectory", "We count tracks starting from 1.");
  if (id > ntracks.size()) Fatal("DrawTrajectory", "There are only %i tracks.", int(ntracks.size()));

  cout << "Start computing tracks..." << endl;

  int counterdaughter = 0;
  bool finished = false;

  int track_id = 0;
  TGeoTrack* curtrack = 0;
  vector<TGeoTrack*> secondary;
  vector<double> X;
  vector<double> Y;
  vector<double> Z;
  double endtime = 0;
  int nframes = 0;

  for (unsigned int i = ntracks[id-1]+1; i < traj.size() && !finished; i++) {
    UHECRTrajectory* tr = traj[i];

    X.push_back(tr->GetPos().X());
    Y.push_back(tr->GetPos().Y());
    Z.push_back(tr->GetPos().Z());

    if (traj[i+1]->GetType() == -1) {
      finished = true;
      endtime = tr->GetTime();
      nframes = i-ntracks[id-1];
    }

    if (int(i) == ntracks[id-1]+1) curtrack = new TGeoTrack(0, int(tr->GetType()));
    if (tr->GetType() == tr->GetInType()) curtrack->AddPoint(tr->GetPos().X(), tr->GetPos().Y(), tr->GetPos().Z(), tr->GetTime());
    else {
      if (secondary.size() && tr->GetType() == secondary.back()->GetPDG()) secondary.back()->AddPoint(tr->GetPos().X(), tr->GetPos().Y(), tr->GetPos().Z(), tr->GetTime());
      else {
	secondary.push_back((TGeoTrack*)curtrack->AddDaughter(counterdaughter, int(tr->GetType())));
	secondary.back()->SetLineColor(counterdaughter);
	secondary.back()->AddPoint(tr->GetPos().X(), tr->GetPos().Y(), tr->GetPos().Z(), tr->GetTime());
	counterdaughter++;
      }
    }
  }
 
    TGeoManager* geom = new TGeoManager("LSS", "The Universe");

  TGeoMaterial *matVacuum = new TGeoMaterial("Vacuum", 0,0,0);
  TGeoMedium *Vacuum = new TGeoMedium("Vacuum 1", 1, matVacuum);

  TGeoVolume *top = geom->MakeBox("TOP", Vacuum, *max_element(X.begin(), X.end()), *max_element(Y.begin(), Y.end()), *max_element(Z.begin(), Z.end()));
  top->SetLineColor(kGreen);
  geom->SetTopVolume(top); // mandatory !

  //--- close the geometry
  geom->CloseGeometry();

  track_id = geom->AddTrack(curtrack);
  geom->SetCurrentTrack(track_id);

  //--- draw the box
  //  geom->SetVisLevel(4);
  top->Draw();

geom->DrawTracks("/*/");
  //geom->AnimateTracks(0,endtime,nframes,"//G/S");

  
  return ;
}

void ReadCRs(char* dirname, std::vector<CRPropaOutput*>& output) {

  string workdir = gSystem->pwd(); // save the current working director

  string format;
  TString filename;
  string fn;
  TList* filelist = NULL;

  TSystemFile sysfile(dirname, NULL); // change to directory where data-files are
  if (sysfile.IsDirectory()) {
    TSystemDirectory sysdir("sysdir", dirname); // change to directory where data-files are
    filelist = sysdir.GetListOfFiles();
    TIter next(filelist);
    
    TObject* obj;
 
    while ( (obj = next()) ) {
      filename = obj->GetName();
      
      if (!filename.CompareTo(".") || !filename.CompareTo("..")) continue;
      if (filename.Contains(".root")) format = "ROOT";
      else if (filename.Contains(".fits")) format = "FITS";
      else format = "ASCII";
      
      fn = std::string(filename.Data());
      
      output.push_back(new CRPropaEvent(fn, format)); // Read in cosmic rays
    }
  }
  else {
    filename = dirname;
    if (!filename.CompareTo(".") || !filename.CompareTo("..")) {
      cerr << "Nothing to do" << endl;
      gSystem->ChangeDirectory(workdir.c_str());
      return ;
    }
    if (filename.Contains(".root")) format = "ROOT";
    else if (filename.Contains(".fits")) format = "FITS";
    else format = "ASCII";
    
    fn = std::string(filename.Data());
    
    output.push_back(new CRPropaEvent(fn, format)); // Read in cosmic rays
  }

  if (filelist) {
    delete filelist;  
    filelist = NULL;
  }
  gSystem->ChangeDirectory(workdir.c_str());
  return ;
}

void ReadNeutrinos(char* dirname, std::vector<CRPropaOutput*>& output) {
  string workdir = gSystem->pwd(); // save the current working director
  string format;
  TString filename;
  string fn;
  TList* filelist = NULL;

  TSystemFile sysfile(dirname, NULL); // change to directory where data-files are
  if (sysfile.IsDirectory()) {
    TSystemDirectory sysdir("sysdir", dirname); // change to directory where data-files are
    filelist = sysdir.GetListOfFiles();
    TIter next(filelist);
    
    TObject* obj;
    
    
    while ( (obj = next()) ) {
      filename = obj->GetName();
      
      if (!filename.CompareTo(".") || !filename.CompareTo("..")) continue;
      if (filename.Contains(".root")) format = "ROOT";
      else if (filename.Contains(".fits")) format = "FITS";
      else format = "ASCII";
      
      fn = std::string(filename.Data());
      
      output.push_back(new CRPropaNeutrino(fn, format)); // Read in neutrinos
    }
  }
  else {
    filename = dirname;
    if (!filename.CompareTo(".") || !filename.CompareTo("..")) {
      cerr << "Nothing to do" << endl;
      gSystem->ChangeDirectory(workdir.c_str());
      return ;
    }
    if (filename.Contains(".root")) format = "ROOT";
    else if (filename.Contains(".fits")) format = "FITS";
    else format = "ASCII";
    
    fn = std::string(filename.Data());
    
    output.push_back(new CRPropaNeutrino(fn, format)); // Read in neutrinos
  }

  if (filelist) {
    delete filelist;  
    filelist = NULL;
  }
  gSystem->ChangeDirectory(workdir.c_str());
  return ;
}

void ReadPhotons(char* dirname, string& options, std::vector<CRPropaOutput*>& output) {
  string workdir = gSystem->pwd(); // save the current working director
  string format;
  TString filename;
  string fn;
  TList* filelist = NULL;

  TSystemFile sysfile(dirname, NULL); // change to directory where data-files are
  if (sysfile.IsDirectory()) {
    TSystemDirectory sysdir("sysdir", dirname); // change to directory where data-files are
    filelist = sysdir.GetListOfFiles();
    TIter next(filelist);
    
    TObject* obj;
    
    
    while ( (obj = next()) ) {
      filename = obj->GetName();
      
      if (!filename.CompareTo(".") || !filename.CompareTo("..")) continue;
      if (filename.Contains(".root")) format = "ROOT";
      else if (filename.Contains(".fits")) format = "FITS";
      else format = "ASCII";
      
      fn = std::string(filename.Data());
      
      output.push_back(new CRPropaPhotonSpectrum(fn, format)); // Read in photons
      if (output.back()->GetGeo() == oneD) options = "HIST P"; // Special option for plotting photons in 1D (switch off statistical errors)
    }
  }
  else {
    filename = dirname;
    if (!filename.CompareTo(".") || !filename.CompareTo("..")) {
      cerr << "Nothing to do" << endl;
      gSystem->ChangeDirectory(workdir.c_str());
      return ;
    }
    if (filename.Contains(".root")) format = "ROOT";
    else if (filename.Contains(".fits")) format = "FITS";
    else format = "ASCII";
    
    fn = std::string(filename.Data());
    
    output.push_back(new CRPropaPhotonSpectrum(fn, format)); // Read in photons
    if (output.back()->GetGeo() == oneD) options = "HIST P"; // Special option for plotting photons in 1D (switch off statistical errors)
  }

  if (filelist) {
    delete filelist;  
    filelist = NULL;
  }
  gSystem->ChangeDirectory(workdir.c_str());
  return ;
}
TH1* ComputeCRSpectrum(const std::vector<CRPropaOutput*>& output, char* histoname, double spect_ind, double old_spect, int nbinsCR, double EminPr, double Emax, double factor) {

  string name;

  TH1D* prot = new TH1D(histoname, ";log_{10}(E / eV);E^{2}dN/dE (eV m^{-2} s^{-1} sr^{-1})", nbinsCR, EminPr, Emax); 
  TH1D* h4 = NULL;

  // compute and add histograms for each run
  for (unsigned int i = 0; i < output.size(); i++) {
    name  = "h4";
    name += i;
    if (output[i]->IsOk()) {
      if (output[i]->WhoAmI() == "CR") {
	h4 = output[i]->ComputeSpectrum(name.c_str(), spect_ind, old_spect, nbinsCR, EminPr, Emax, factor);
	prot->Add(h4);
      }
      if (h4) delete h4;
      h4 = NULL;
    }
  }
  
  return prot;
}

TH1* ComputeMassSpectrum(const std::vector<CRPropaOutput*>& output, char* histoname, double spect_ind, double old_spect, int nbinsCR, double EminPr, double Emax) {

  string msname;

  TProfile* ms = new TProfile(histoname, ";log_{10}(E / eV);< A >", nbinsCR, EminPr, Emax); 

  TProfile* h3 = NULL;

  // compute and add histograms for each run
  for (unsigned int i = 0; i < output.size(); i++) {
    if (output[i]->IsOk()) {
      if (output[i]->WhoAmI() == "CR") {
	msname = "h3";
	msname += i;
	h3 = output[i]->ComputeMassSpectrum(msname.c_str(), spect_ind, old_spect, nbinsCR, EminPr, Emax);
	ms->Add(h3);
	if (h3) {
	  delete h3;
	  h3 = NULL;
	}
      }
    }
  }
  
  return ms;
}

TH1* ComputeNeutrinoSpectrum(const std::vector<CRPropaOutput*>& output, char* histoname, double spect_ind, double old_spect, int nbins, double Emin, double Emax, double factor, bool& found_neut) {
  string name;

  TH1D* neutrinos = new TH1D(histoname, ";log_{10}(E / eV);E^{2}dN/dE (eV m^{-2} s^{-1} sr^{-1})", nbins, Emin, Emax);
  neutrinos->SetMarkerColor(kBlue);

  TH1D* h4 = NULL;

  // compute and add histograms for each run
  for (unsigned int i = 0; i < output.size(); i++) {
    name  = "h4";
    name += i;
    if (output[i]->IsOk()) {
      if (output[i]->WhoAmI() == "Neutrino") { 
	h4 = output[i]->ComputeSpectrum(name.c_str(), spect_ind, old_spect, nbins, Emin, Emax, factor);
	neutrinos->Add(h4); 
	found_neut = true; 
	if (h4) {
	  delete h4;
	  h4 = NULL;
	}
      }
    }
  }

  return neutrinos;
}

TH1* ComputePhotonSpectrum(const std::vector<CRPropaOutput*>& output, char* histoname, double spect_ind, double old_spect, int nbins, double Emin, double Emax, double factor, double reweightfactor, bool& found_phot) {
  string name;

  TH1D* photons = new TH1D(histoname, ";log_{10}(E / eV);E^{2}dN/dE (eV m^{-2} s^{-1} sr^{-1})", nbins, Emin, Emax); 
  
  TH1D* h4 = NULL;

  // compute and add histograms for each run
  for (unsigned int i = 0; i < output.size(); i++) {
    name  = "h4";
    name += i;
    if (output[i]->IsOk()) {
      if (output[i]->WhoAmI() == "PhotonSpectrum") { 
	h4 = output[i]->ComputeSpectrum(name.c_str(), spect_ind, old_spect, nbins, Emin, Emax, factor, reweightfactor);
	photons->Add(h4); 
	found_phot = true; 
	if (h4) {
	  delete h4;
	  h4 = NULL;
	}
      }
    }
  }

  return photons;
}

TProfile* ComputeDeflections(const std::vector<CRPropaOutput*>& output, const char* histoname, const double& spect_ind, const double& old_spect, const int& nbins, const double& Emin, const double& Emax) {

  string name;
  TProfile* deflections = new TProfile(histoname, ";log_{10}(E / eV);#delta (#circ)", nbins, Emin, Emax); 

  TProfile* h3 = NULL;

  for (unsigned int i = 0; i < output.size(); i++) {
    if (output[i]->IsOk()) {
      if (output[i]->WhoAmI() == "CR") {
	h3 = output[i]->ComputeDeflections(name.c_str(), spect_ind, old_spect, nbins, Emin, Emax);
	deflections->Add(h3);
	h3->Delete();

      }
    }
  }


  return deflections;
}

TProfile* ComputeCumulativeDeflections(const std::vector<CRPropaOutput*>& output, const char* histoname, const double& spect_ind, const double& old_spect, const int& nbins, const double& Emin, const double& Emax, const double& Dmax) {

  string name;

  TProfile* cumdefl = new TProfile(histoname, ";#delta (#circ);", nbins, 0, 180); 

  TProfile* h3 =NULL;
  for (unsigned int i = 0; i < output.size(); i++) {
    if (output[i]->IsOk()) {
      if (output[i]->WhoAmI() == "CR") {
	h3 = output[i]->ComputeCumulativeDeflections(name.c_str(), spect_ind, old_spect, nbins, Emin, Emax, Dmax);
	cumdefl->Add(h3);
	h3->Delete();

      }
    }
  }

  return cumdefl;
}

float* ComputeSkyMap(const std::vector<CRPropaOutput*>& output, const double& spect_ind, const double& old_spect, const double& Emin, const double& Emax, const double& Dmax, const long& nside) {

  unsigned long npix = 12*nside*nside;

  float* skymap = new float[npix];

  for (unsigned int i = 0; i < output.size(); i++) {
    if (output[i]->IsOk()) {
      if (output[i]->WhoAmI() == "CR") {
	vector<double> fakemap = output[i]->ComputeSkyMap(spect_ind, old_spect, Emin, Emax, Dmax, nside);
	for (unsigned long j = 0; j < npix; ++j) {
	  skymap[j] += fakemap[j];
	}
      }
    }
  }

  return skymap;

}
