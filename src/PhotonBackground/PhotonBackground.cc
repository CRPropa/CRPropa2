#include <PhotonBackground.h>

TGraph* RootPlotSpectrum(PhotonBackground& BGrd,
			 double Emin,
			 double Emax,
			 int EBins,
			 double x,
			 double y,
			 double z,
			 double redshift){
  
  std::cout<<"Emin="<<Emin<<"\tEmax="<<Emax<<std::endl;

  std::vector<double> E_Vec;
  std::vector<double> PhotD_Vec;
  
  for(int i=0; i<EBins; i++){
    double E = log10(Emin)+i*fabs(log10(Emax)-log10(Emin))/(double) EBins;
    double PhotonDensity = BGrd.GetPhotonDensity(x,
						 y,
						 z,
						 redshift,
						 pow(10.,E));
    std::cout<<"E="<<E<<"\tPhotD="<<PhotonDensity<<std::endl;
    E_Vec.push_back(E);
    PhotD_Vec.push_back(PhotonDensity);
  } 
  
  TGraph* TheGraph = new TGraph(E_Vec.size(), &(E_Vec[0]), &(PhotD_Vec[0]));
  
  return(TheGraph);
}
 
