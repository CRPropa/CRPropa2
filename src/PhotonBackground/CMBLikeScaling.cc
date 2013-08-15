#include <CMBLikeScaling.h> 

CMBLikeScaling:: CMBLikeScaling(std::string XMLFileName): TIRBzEvolutionModel::TIRBzEvolutionModel(XMLFileName){
  //Read the MaxRedshift value from xml configuration
  bool lCheckIR = XmlExtract().CheckElement("InfraredBackground") ;
  if (lCheckIR) {
    TiXmlElement* lpXmlIR = XmlExtract().GetElement("InfraredBackground") ;
    if (lpXmlIR) {
      TiXmlElement* lpXmlZmax = lpXmlIR->FirstChildElement("MaxRedshift") ;
      if(lpXmlZmax) lpXmlZmax->Attribute("value",&_fMaxRedshift) ;
    } else throw TXmlErr("No MaxRedshift given.") ;
  }   
  std::cout<<"MaxRedshift for CMB-like scaling is: "<<_fMaxRedshift<<std::endl;
}

double CMBLikeScaling::GetScalingFactor(double redshift) {
  if(redshift>_fMaxRedshift) return(0.);
 
  return pow(1+redshift,3);
}
