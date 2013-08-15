#include "IRBzEvolutionFactory.h"
#include <new>

TIRBzEvolutionModel* IRBzEvolutionFactory::_fModelInstance=NULL;
IRBzEvolutionFactory* IRBzEvolutionFactory::_fFactoryInstance=NULL;
std::string IRBzEvolutionFactory::_fFileName="";

/////////////////////////////////////////////////////////////////////////
//FACTORY CLASS
//Method to create and call the Factory after it was already constructed unsing the name of the xml file.
IRBzEvolutionFactory* IRBzEvolutionFactory::GetInstance(std::string XMLFileName){
  if(! _fFactoryInstance){
    _fFactoryInstance    = new IRBzEvolutionFactory (XMLFileName);
    _fModelInstance	 = IRBzEvolutionFactory::GetConcreteIRBzEvolutionModel(XMLFileName);
    _fFileName           = XMLFileName;
  }
  //We deal with a singleton -> configuration is not allowed to be changed
  if(XMLFileName.compare(_fFileName)!=0) throw TXmlErr("IRBzEvolutionFactory::GetInstance(...): path to xml file has changed.");
      return _fFactoryInstance;
}
//If the singleton factory has already been constructed using the xml file location (s. GetInstance(XMLFileName) above)  this method without the xml-file name can be used, too. This is needed as e.g. TNucleus doe not know the xml configuration path but needs the IRB scaling in spite of that.  
IRBzEvolutionFactory* IRBzEvolutionFactory::GetInstance(){
  if(! _fFactoryInstance){
    throw TXmlErr("IRBzEvolutionFactory::GetInstance() called but the first call need to give the path to the xml file.");
  }
  return _fFactoryInstance;
}

TIRBzEvolutionModel* IRBzEvolutionFactory::GetConcreteIRBzEvolutionModel(std::string XMLFileName){

  //Read settings from xml configuration
  //TXmlParam Test(XMLFileName.c_str());
  bool lCheckIRBzModel = IRBzEvolutionFactory::GetInstance()->XmlExtract().CheckElement("IRB_MFPScalingModel") ;
  if ( ! lCheckIRBzModel){
    std::cout<<"No IRB evolution scaling model chosen via xml settings. Default \"cmblike\"-scaling will be used."<<std::endl;
#ifdef UBUNTU 
    return(new typename CMBLikeScaling::CMBLikeScaling(XMLFileName));
#else 
    return(new CMBLikeScaling(XMLFileName));
#endif
  }else{
    TiXmlElement* lpXmlIRBzModel = IRBzEvolutionFactory::GetInstance()->XmlExtract().GetElement("IRB_MFPScalingModel") ;
    std::string lIRBzModelType = lpXmlIRBzModel->Attribute("type") ;
    if (lIRBzModelType == "CMBLike") {
      std::cout<<"IRB_MFPScalingModel: "<<lIRBzModelType<<"is selected."<<std::endl;
#ifdef UBUNTU
      return(new typename CMBLikeScaling::CMBLikeScaling(XMLFileName));
#else
      return(new CMBLikeScaling(XMLFileName));
#endif
    }else if (lIRBzModelType == "Kneiske2004_BestFit") {
      std::cout<<"IRB_MFPScalingModel: "<<lIRBzModelType<<" is selected."<<std::endl;
#ifdef UBUNTU
      return(new typename Kneiske2004_BestFit::Kneiske2004_BestFit(XMLFileName));
#else
      return(new Kneiske2004_BestFit(XMLFileName));
#endif
    }else if (lIRBzModelType == "UserDefinedScaling") {
	std::cout<<"IRB_MFPScalingModel: "<<lIRBzModelType<<" is selected."<<std::endl;
#ifdef UBUNTU
	return(new typename UserDefinedScaling::UserDefinedScaling(XMLFileName));
#else
	return(new UserDefinedScaling(XMLFileName));
#endif
  }else if (lIRBzModelType == "None") {
  std::cout<<"IRB_MFPScalingModel: "<<lIRBzModelType<<" is selected."<<std::endl;
#ifdef UBUNTU
      return(new typename NoIRBScaling::NoIRBScaling(XMLFileName));
#else
      return(new NoIRBScaling(XMLFileName));
#endif
    }else  throw TCrpErr("<IRB_MFPScalingModel type=XYZ> selected but type not specified.") ;
  }
}

void IRBzEvolutionFactory::PlotAllModels(){
  double zMin=0., zMax=7.;
  int NBins=10000;
  
  TGraph Kneiske2004_BestFit_TG      =  Kneiske2004_BestFit(_fFileName).PlotModelTGraph(zMin, zMax, NBins);
  TGraph CMBLikeScaling_TG   =  CMBLikeScaling(_fFileName).PlotModelTGraph(zMin, zMax, NBins);
  TGraph NoIRBScaling_TG     =  NoIRBScaling(_fFileName).PlotModelTGraph(zMin, zMax, NBins);

  TCanvas AllIRBzEvolutionModelsCan("AllIRBzEvolutionModelsCan"
				    ,"AllIRBzEvolutionModelsCan"
				    ,1);
  TLegend legend(.75,.80,.95,.95);
  int a=0;
  Kneiske2004_BestFit_TG.Draw("Al");
  legend.AddEntry(&Kneiske2004_BestFit_TG,"Kneiske 2004 (best fit model)");
  Kneiske2004_BestFit_TG.SetLineStyle(++a);
  CMBLikeScaling_TG.Draw("l");
  legend.AddEntry(&CMBLikeScaling_TG,"cmb-like scaling + cut-off z>z_{max}");
  CMBLikeScaling_TG.SetLineStyle(++a);
  NoIRBScaling_TG.Draw("l");
  legend.AddEntry(&NoIRBScaling_TG,"no scaling");
  NoIRBScaling_TG.SetLineStyle(++a);
  legend.Draw("same");
  AllIRBzEvolutionModelsCan.SaveAs("AllIRBzEvolutionModelsCan.root");
}
