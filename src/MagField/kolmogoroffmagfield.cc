/**
   @file    kolmogoroffmagfield.cc
   @author  Eric Armengaud, armengau@in2p3.fr
   @brief   Implementation of the TKolmogoroffMagfield class. See the .h file
*/

#include "kolmogoroffmagfield.h"
#include "fftw3.h"

TKolmogoroffMagField::TKolmogoroffMagField(const char* aFileName) : TXmlParam(aFileName) {
  SetType(MAGFIELD_KOLMOGOROFF);

  // xml initialization
  string lDum ;
  TiXmlElement* lpXmlField = XmlExtract().GetElement("MagneticField") ;
  TiXmlElement* lpXmlNx = lpXmlField->FirstChildElement("Nx") ;
  TiXmlElement* lpXmlNy = lpXmlField->FirstChildElement("Ny") ;
  TiXmlElement* lpXmlNz = lpXmlField->FirstChildElement("Nz") ;
  TiXmlElement* lpXmlStep = lpXmlField->FirstChildElement("Step_Mpc") ;
  if (!lpXmlNx || !lpXmlNy || !lpXmlNz || !lpXmlStep)
    throw TCrpErr("Magnetic grid not properly specified in xml file");
  lDum = lpXmlNx->Attribute("value",&_fNx) ;
  lDum = lpXmlNy->Attribute("value",&_fNy) ;
  lDum = lpXmlNz->Attribute("value",&_fNz) ;
  lDum = lpXmlStep->Attribute("value",&_fStepsize) ;
  _fStepsize *= Mpc ;
  TiXmlElement* lpXmlOrigin = lpXmlField->FirstChildElement("Origin") ;
  if (lpXmlOrigin) {
    TiXmlElement* lpXmlOx = lpXmlOrigin->FirstChildElement("X_Mpc") ;
    TiXmlElement* lpXmlOy = lpXmlOrigin->FirstChildElement("Y_Mpc") ;
    TiXmlElement* lpXmlOz = lpXmlOrigin->FirstChildElement("Z_Mpc") ;
    if (!lpXmlOx || !lpXmlOy || !lpXmlOz)
      throw TCrpErr("Magnetic grid origin not properly specified");
    double lOx,lOy,lOz;
    lDum = lpXmlOx->Attribute("value",&lOx) ;
    lDum = lpXmlOy->Attribute("value",&lOy) ;
    lDum = lpXmlOz->Attribute("value",&lOz) ;
    _fOrigin.set(lOx*Mpc,lOy*Mpc,lOz*Mpc) ;
  } else _fOrigin.set(0,0,0) ;


  TiXmlElement* lpXmlSpec = lpXmlField->FirstChildElement("SpectralIndex") ;
  if (!lpXmlSpec) throw TCrpErr("Kolmogoroff B-field requires spectral index.") ;
  lDum = lpXmlSpec->Attribute("value",&_fSpectralIndex) ;

  TiXmlElement* lpXmlKmin = lpXmlField->FirstChildElement("Kmin") ;
  TiXmlElement* lpXmlKmax = lpXmlField->FirstChildElement("Kmax") ;
  if (!lpXmlKmin || !lpXmlKmax) throw TCrpErr("Kolmogoroff B-field requires Kmax and Kmin.") ;
  else {
    lDum = lpXmlKmin->Attribute("value",&_fKmin) ;
    lDum = lpXmlKmax->Attribute("value",&_fKmax) ;
  }

  if (_fNx != _fNy || _fNy != _fNz) throw TCrpErr("Kolmogoroff B-field requires Nx=Ny=Nz.") ;
  unsigned int N = _fNx ;
  unsigned int N2 = int(_fNx)/2+1 ;


  std::cout<<"Kolmogorov Field"<<std::endl;
  std::cout<<"Kmin "<<_fKmin<<std::endl;
  std::cout<<"Kmax "<<_fKmax<<std::endl;

  double lBrms ; // sqrt( <B(k)^2> ) = sqrt( <B(x)^2> )
  TiXmlElement* lpXmlAmpl = lpXmlField->FirstChildElement("RMS_muG") ;
  lDum = lpXmlAmpl->Attribute("value",&lBrms) ;
  if (!lpXmlAmpl) throw TCrpErr("B field RMS not specified.");
  lBrms *= muG ;
  // end of xml intialization


  // create field
  // vector components of B-field
  fftw_complex *Bx, *By, *Bz;
  Bx = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N*N*N);
  By = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N*N*N);
  Bz = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N*N*N);


  // possible wavenumbers for a grid of N steps
  float K[N];
  for (int i=0; i<N/2; i++) {
    K[i] = (float)i/N;
  }
  for (int i=N/2; i<N; i++) {
    K[i] = (float)i/N-1;
  }

  int i;
  double k, theta, phase, cosPhase, sinPhase;

  TVector3D e1, e2, ek, b;
  TVector3D n0(1,1,1); // arbitrary vector for orthogonal construction

  for (int ix=0; ix<N; ix++) {
    for (int iy=0; iy<N; iy++) {
      for (int iz=0; iz<N; iz++) {

        int i = ix*N*N + iy*N + iz;
        ek.set(K[ix], K[iy], K[iz]);
        k = ek.mag();
        if ((k < _fKmin) || (k > _fKmax)) {
			Bx[i][0] = 0;
			Bx[i][1] = 0;
			By[i][0] = 0;
			By[i][1] = 0;
			Bz[i][0] = 0;
			Bz[i][1] = 0;
			continue;
        		}
        else // wave inside turbulent range
        {
          // construct an orthogonal base e1, e2, ek
          if ((ix == iy) && (iy == iz)) // ek || (1,1,1)
          {
            e1.set(-1.,1.,0.);
            e2.set(1.,1.,-2.);
          }
          else // ek not || (1,1,1)
          {
            e1 = n0.cross(ek);
            e2 = ek.cross(e1);
          }
          e1 /= e1.mag();
          e2 /= e2.mag();

          // random orientation perpendicular to k
          theta = 2 * M_PI *RandFlat::shoot();
          b = e1*cos(theta) + e2*sin(theta);

          // gaussian amplitude weighted with k^alpha/2
          b *= RandGauss::shoot() * pow(k, _fSpectralIndex/2.);

          // uniform random phase
          phase = 2 * M_PI *RandFlat::shoot();
          fftw_complex a[2];
          cosPhase = cos(phase); // real part
          sinPhase = sin(phase); // imaginary part

          Bx[i][0] = b.x() * cosPhase; Bx[i][1] = b.x() * sinPhase;
          By[i][0] = b.y() * cosPhase; By[i][1] = b.y() * sinPhase;
          Bz[i][0] = b.z() * cosPhase; Bz[i][1] = b.z() * sinPhase;
        }
      }
    }
  }

  // perform inverse FFT on each component
  fftw_plan p;
  p = fftw_plan_dft_3d(N, N, N, Bx, Bx, FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_execute(p);
  p = fftw_plan_dft_3d(N, N, N, By, By, FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_execute(p);
  p = fftw_plan_dft_3d(N, N, N, Bz, Bz, FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_execute(p);
  fftw_destroy_plan(p);

  // normalize RMS to 1
  double norm = 0;
  for (unsigned int ix=0; ix<N; ix++)
    for (unsigned int iy=0; iy<N; iy++)
      for (unsigned int iz=0; iz<N; iz++)
      {
        i = ix*N*N + iy*N + iz;
        norm +=  pow(Bx[i][0],2.)+pow(By[i][0],2.)+pow(Bz[i][0],2.) ;
      }
  norm = sqrt(norm/(N*N*N));
  for (unsigned int ix=0; ix<N; ix++)
    for (unsigned int iy=0; iy<N; iy++)
      for (unsigned int iz=0; iz<N; iz++)
      {
        i = ix*N*N + iy*N + iz;
        Bx[i][0] /= norm;
        By[i][0] /= norm;
        Bz[i][0] /= norm;
      }

//  // output
//  for (int ix=0; ix<N; ix++)
//    for (int iy=0; iy<N; iy++)
//      for (int iz=0; iz<N; iz++)
//      {
//        i = ix*N*N + iy*N + iz;
//        cout << Bx[i][0] << "," << By[i][0] << "," << Bz[i][0] << endl;
//      }

  // make field usable for crpropa
  _fpB = new TVector3D[N*N*N];
  for (unsigned int ix=0; ix<N; ix++)
    for (unsigned int iy=0; iy<N; iy++)
      for (unsigned int iz=0; iz<N; iz++)
        {
        i = ix*N*N + iy*N + iz;
        _fpB[i].setX(Bx[i][0] * lBrms);
        _fpB[i].setY(By[i][0] * lBrms);
        _fpB[i].setZ(Bz[i][0] * lBrms);


	//        cout<<_fpB[i]<<endl;
        }

  _fBmax = 0 ;
  for (int i=0; i<_fNx; i++)
    for (int j=0; j<_fNy; j++)
      for (int k=0; k<_fNz; k++)
        _fBmax = max(_fBmax, (_fpB[i*_fNy*_fNz+j*_fNz+k]).mag() ) ;

  fftw_free(Bx); fftw_free(By); fftw_free(Bz);

}


TKolmogoroffMagField::~TKolmogoroffMagField() {
  delete[] _fpB ;
}
