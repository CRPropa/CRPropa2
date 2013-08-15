#include "uniforminfrared.h"


/* Remember to handle unit conversions !!!!! */

const double FLUX_CONVERSION = 3.82182e3;

TUniformInfrared::TUniformInfrared() : IR() {
    SetType(IR_UNIFORM_PRIMACK);
    //PRIMACK
    std::vector<double> XData;
    XData.push_back(-1.); 
    XData.push_back(-0.75); 
    XData.push_back(-0.5); 
    XData.push_back(-0.25); 
    XData.push_back(0.);
    XData.push_back(0.25); 
    XData.push_back(0.5); 
    XData.push_back(0.75); 
    XData.push_back(1.); 
    XData.push_back(1.25); 
    XData.push_back(1.5); 
    XData.push_back(1.75); 
    XData.push_back(2.); 
    XData.push_back(2.25); 
    XData.push_back(2.5);
    
    std::vector<double> YData;
    YData.push_back(0.8); 
    YData.push_back(1.1); 
    YData.push_back(1.15); 
    YData.push_back(1.2); 
    YData.push_back(1.3);
    YData.push_back(1.2); 
    YData.push_back(1.05); 
    YData.push_back(0.7); 
    YData.push_back(0.4); 
    YData.push_back(0.3); 
    YData.push_back(0.5); 
    YData.push_back(0.8); 
    YData.push_back(1.1); 
    YData.push_back(1.3); 
    YData.push_back(1.);
    
    spectrum = new TPhotonSpectrum(XData,YData);
}

TUniformInfrared::TUniformInfrared(const char* aFileName) : IR(), TXmlParam(aFileName) {

  string lDum ;
  string lSpectralShape = "Primack" ;
  //  _fZmax = DEFAULT_ZMAX_IR ;

  bool lCheckIR = XmlExtract().CheckElement("InfraredBackground") ;
  if (lCheckIR) {
    TiXmlElement* lpXmlIR = XmlExtract().GetElement("InfraredBackground") ;
    if (lpXmlIR) {
      TiXmlElement* lpXmlSpectrum = lpXmlIR->FirstChildElement("SpectralShape") ;
      if (lpXmlSpectrum) lSpectralShape = lpXmlSpectrum->Attribute("type") ;
      TiXmlElement* lpXmlZmax = lpXmlIR->FirstChildElement("MaxRedshift") ;
      if (lpXmlZmax) lDum = lpXmlZmax->Attribute("value",&_fZmax) ; 
    }
  }

  if (lSpectralShape == "Primack" || lSpectralShape == "Kneiske") {
    //The 
    SetType(IR_UNIFORM_PRIMACK);
    //PRIMACK (The spectra is actually Primack, but it is nowhere used...) 
    std::vector<double> XData;
    XData.push_back(-1.); 
    XData.push_back(-0.75); 
    XData.push_back(-0.5); 
    XData.push_back(-0.25); 
    XData.push_back(0.);
    XData.push_back(0.25); 
    XData.push_back(0.5); 
    XData.push_back(0.75); 
    XData.push_back(1.); 
    XData.push_back(1.25); 
    XData.push_back(1.5); 
    XData.push_back(1.75); 
    XData.push_back(2.); 
    XData.push_back(2.25); 
    XData.push_back(2.5);
    
    std::vector<double> YData;
    YData.push_back(0.8); 
    YData.push_back(1.1); 
    YData.push_back(1.15); 
    YData.push_back(1.2); 
    YData.push_back(1.3);
    YData.push_back(1.2); 
    YData.push_back(1.05); 
    YData.push_back(0.7); 
    YData.push_back(0.4); 
    YData.push_back(0.3); 
    YData.push_back(0.5); 
    YData.push_back(0.8); 
    YData.push_back(1.1); 
    YData.push_back(1.3); 
    YData.push_back(1.);
    
    spectrum = new TPhotonSpectrum(XData,YData);
    
  } else if (lSpectralShape == "High") {
    SetType(IR_UNIFORM_HIGH);
  } else if (lSpectralShape == "Low") {
    SetType(IR_UNIFORM_LOW);
  } else throw TXmlErr("Unknown IR background spectral type") ;


}

double TUniformInfrared::GetPhotonDensity(double x, 
					  double y, 
					  double z, 
					  double redshift,
					  double PhotonEnergy) {
  
  double REDSHIFT, Z, ZMAX_IR=1.;
  //Redshift evolution as for CMB. (added Dec.'05)
  //conversion from nW/cm^3/sr to eV/cm^3
  //        FLUX_CONVERSION = 2.9979e10/4./dacos(-1.)*1.602e-19*1.e9*1.e4
  
  //c conversion from eV to micrometers
  double X = 1.2398*(1.+redshift)/PhotonEnergy;
  double RESULT=0.;

  vector<double> XData(spectrum->GetE());
  vector<double> YData(spectrum->Spectrum());


  /* Check with others */
  if (X>0.0){
    RESULT=0.;
  }   
  else if(log10(X)>=XData[15-1]){ 
    RESULT = (YData[15-1] - YData[14-1])/(XData[15-1] - XData[14-1])*(log10(X) - XData[14-1]) + YData[14-1];
    RESULT = pow(10.,RESULT);
  }
  
  if (log10(X) <= XData[1-1]) RESULT=0.;
  
  if ((log10(X) < XData[15-1]) && (log10(X) > XData[1-1])){
    int INDEX = 2;
    do 
      INDEX = INDEX+1;
    while (XData[INDEX-1] < log10(X));
    
    RESULT = (YData[INDEX-1] - YData[INDEX-1-1])/ (XData[INDEX-1] - XData[INDEX-1-1])*(log10(X) - XData[INDEX-1-1])+  YData[INDEX-1-1];
    
    RESULT = pow(10., RESULT);  
  }
  
  if (Z>ZMAX_IR) {
    return(0.);
  }
  
#ifdef DEBUG_OUTPUT
  //std::cout<<"Density="
  //   <<RESULT*pow((1.+redshift), 2)/pow((PhotonEnergy/(1.+redshift)), 2)/FLUX_CONVERSION
  //   <<std::endl;
#endif
  return(RESULT*pow((1.+redshift), 2)/pow((PhotonEnergy/(1.+redshift)), 2)/FLUX_CONVERSION);



}

//TPhotonSpectrum TUniformInfrared::Spectrum(double x, double y = -1, double z = -1) { return *spectrum;}

