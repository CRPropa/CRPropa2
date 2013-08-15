//////////////////////////////////////////////////////////////////////////////////
//ABSTRACT MODEL CLASS (from which the concrete IRB redshift evolution models are derived).
#include <TIRBzEvolutionModel.h>

TIRBzEvolutionModel::TIRBzEvolutionModel(std::string XMLFileName): TXmlParam(XMLFileName.c_str()), _fMaxRedshift(DEFAULT_ZMAX_IR){}

TGraph TIRBzEvolutionModel::PlotModelTGraph(double zMin, double zMax, int NPoints){
  
  std::vector<double> x,y;
  for(double z=zMin; z<=zMax; z+=(zMax-zMin)/( (double) NPoints )){
    x.push_back(z);
    y.push_back(this->GetScalingFactor(z));
  }
  
  return( TGraph( x.size(), &(x[0]), &(y[0]) ) );
}
