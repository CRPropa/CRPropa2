
/**
   @file    gridfield.cc
   @author  Eric Armengaud, armengau@in2p3.fr
   @brief   Implementation of the TGridField class. See the .h file
*/

#include "gridfield.h"

TVector3D TGridField::getField(TVector3D aPosition) const {

  // aPosition must be inside a shell around the grid
  TVector3D lLocB ;
  int lQx,lQy,lQz ;
  int lQx1,lQy1,lQz1 ;
  double a,b,c,d,e,f ;
  TVector3D lPosInStep = (aPosition-_fOrigin)/_fStepsize ;
  lQx = (int)floor(lPosInStep.x()) ;
  a=lPosInStep.x()-lQx ;
  d=1-a ;
  lQx1 = lQx+1 ;
  lQx = lQx % _fNx ;
  if (lQx<0) lQx += _fNx ;
  lQx1 = lQx1 % _fNx ;
  if (lQx1<0) lQx1 += _fNx ;
  lQy = (int)floor(lPosInStep.y()) ;
  b=lPosInStep.y()-lQy ;
  e=1-b ;
  lQy1 = lQy+1 ;
  lQy = lQy % _fNy ;
  if (lQy<0) lQy += _fNy ;
  lQy1 = lQy1 % _fNy ;
  if (lQy1<0) lQy1 += _fNy ;
  lQz = (int)floor(lPosInStep.z()) ;
  c=lPosInStep.z()-lQz ;
  f=1-c ;
  lQz1 = lQz+1 ;
  lQz = lQz % _fNz ;
  if (lQz<0) lQz += _fNz ;
  lQz1 = lQz1 % _fNz ;
  if (lQz1<0) lQz1 += _fNz ;
  
  int NyNz = _fNy*_fNz; 
  lLocB=d*
    (e*(f*_fpB[lQx*NyNz+lQy*_fNz+lQz]+c*_fpB[lQx*NyNz+lQy*_fNz+lQz1])+
     b*(f*_fpB[lQx*NyNz+lQy1*_fNz+lQz]+c*_fpB[lQx*NyNz+lQy1*_fNz+lQz1]))+
    a*(e*(f*_fpB[lQx1*NyNz+lQy*_fNz+lQz]+c*_fpB[lQx1*NyNz+lQy*_fNz+lQz1])+
       b*(f*_fpB[lQx1*NyNz+lQy1*_fNz+lQz]+c*_fpB[lQx1*NyNz+lQy1*_fNz+lQz1])) ;
  
  if ( lLocB.mag() > this->Bmax() ) 
    throw TCrpErr("Interpolated B field larger than Bmax.") ;
  
  return lLocB;
}


int TGridField::writeMagField(string outFName){


      fitsfile *foutFile;
      int status=0; 
      char errorstr[10];
      
      fits_create_file(&foutFile, outFName.c_str(), &status); 
      if(status!=0){
	fits_report_error(stderr, status);
	sprintf(errorstr, "%i", status);
	throw TCrpErr( string("Error creating magnetic field output file.") + string(errorstr));
      }

      static const int bitpix=-32;
      static const int naxis=3;
      long naxes=256;

       long Nx=this->Nx();
       long Ny=this->Nx();
       long Nz=this->Nz();
       long Nvec[3]={Nx, Ny, Nz};// Convenience 
      TVector3D* magField;
      magField=this->B();
      double *Bx=new double[Nx*Ny*Nz];
      double *By=new double[Nx*Ny*Nz];
      double *Bz=new double[Nx*Ny*Nz];


      unsigned long MaxInd=Nx*Ny*Nz;
      unsigned long lInd, lIndFits;
      for(unsigned long i=0; i< Nx; i++){
	for(unsigned long j=0; j< Ny; j++){
	  for(unsigned long k=0; k<Nz; k++){
	    lIndFits=i+j*Nx+k*Nx*Ny;
	    lInd= i*Ny*Nz+ j*Nz+k;
	    
	    Bx[lIndFits]=(magField+lInd)->x()/gauss;
	    By[lIndFits]=(magField+lInd)->y()/gauss;
	    Bz[lIndFits]=(magField+lInd)->z()/gauss;
	    
	  }
	}
      }

      int HDUnum=-3;
      fits_get_hdu_num(foutFile, &HDUnum);
      int firstPixel=1;

      fits_create_img(foutFile, FLOAT_IMG, 3, Nvec, &status);
      if(status!=0){
	std::cout<< "Error in Bx header" << std::endl;
	fits_report_error(stderr, status);
	exit(status);
      }
      fits_write_img(foutFile, TDOUBLE, firstPixel, Nx*Ny*Nz, Bx, &status);
      if(status!=0){
	std::cout<< "Error in Bx writing" << std::endl;
	fits_report_error(stderr, status);
	exit(status);
      }
      fits_create_img(foutFile, FLOAT_IMG, 3, Nvec, &status);
      fits_write_img(foutFile, TDOUBLE, firstPixel, Nx*Ny*Nz, By, &status);
      if(status!=0){
	std::cout<< "Error in By" << std::endl;
	fits_report_error(stderr, status);
	exit(status);
      }
      fits_create_img(foutFile, FLOAT_IMG, 3, Nvec, &status);
      fits_write_img(foutFile, TDOUBLE, firstPixel, Nx*Ny*Nz, Bz, &status);
      if(status!=0){
	std::cout<< "Error in Bx" << std::endl;
	fits_report_error(stderr, status);
	exit(status);
      }


      fits_close_file(foutFile, &status);
      if(status!=0){
	fits_report_error(stderr, status);
	exit(status);
      }

      return 1;

}
