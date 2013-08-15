
/**
   @file    everything.cc
   @author  Tristan Beau, beau@in2p3.fr
   @brief   Implementation Class which will contain every usefull info for simu
*/

/*changed by Nils Nierstenhoefer (NN) to add Interactions type=nucleus*/

#include "everything.h"

TEveryThing::TEveryThing(const char *aFileName ) {

  Basic = new TBasicParam( aFileName );
  TXmlExtract lExtract(aFileName) ;
  TiXmlElement* lpXmlUniverse = lExtract.GetElement("Environment") ;
  string lUnivType = lpXmlUniverse->Attribute("type") ;
  if (lUnivType == "One Dimension") {
    Univ = new TEnv1D(aFileName) ;
  }  else if (lUnivType == "LSS")
    Univ = new TLargeScaleStructure(aFileName) ;
  else if (lUnivType == "Galactic") 
    Univ = new TGalacticStructure(aFileName) ;
  else
    throw TXmlErr("Incorrect environment type.");

  // Consistency tests
  if ( Basic->RecordMode() == "Full Trajectories" ) {
    if ( Univ->InteractionData()->SecPhotonFlag() )
      throw TXmlErr("No photon secondaries in Full Trajectories mode.") ;
    if ( Univ->InteractionData()->SecNuFlag() )
      throw TXmlErr("No neutrino secondaries in Full Trajectories mode.") ;
  }
 
  if ( Basic->Emin() < 0.07*EeV && Univ->InteractionData()->PairProdFlag() 
       && !(Univ->Sources()->IsPhotonSources()) ) // Photon sources : allow small Emin
    //throw TXmlErr("Emin < 0.07 EeV : incompatible with pair production.") ;
    std::cout<<" Warning: Emin < 0.07 EeV : incompatible with pair production. "<<std::endl;
  if ( Basic->Emin() < 10*EeV && Univ->InteractionData()->Type() == INTERACTION_BASIC 
       && !(Univ->Sources()->IsPhotonSources()) ) // Photon sources : allow small Emin
    throw TXmlErr("Emin < 10 EeV : incompatible with F77 interactions.") ;
  if ( Univ->Sources()->IsPhotonSources() && Basic->Emin() < 10.*MeV )
    throw TXmlErr("Emin < 10 MeV for photons : incompatible with DINT.") ;

  if ( Univ->Sources()->Ecut() >= 1000*ZeV && Univ->InteractionData()->PairProdFlag() )
    throw TXmlErr("Emax > 1000 ZeV : incompatible with pair production.") ;
  if ( Univ->Sources()->Ecut() >= 100*ZeV && Univ->InteractionData()->Type() == INTERACTION_BASIC )
    throw TXmlErr("Emax > 1000 ZeV : incompatible with F77 interactions.") ;

  if ( Univ->Sources()->IsPhotonSources() && Univ->Type() != UNIVERSE_ENV1D )
    throw TXmlErr("No direct photon injection is possible if the environment is not 1D.") ;
  if ( Univ->Sources()->IsPhotonSources() && Univ->InteractionData()->Type() != INTERACTION_PHOTON)
    throw TXmlErr("Interactions should be of photon type with photon sources.") ;
  if ( !(Univ->Sources()->IsPhotonSources()) && 
       Univ->InteractionData()->Type() == INTERACTION_PHOTON)
    throw TXmlErr("No single photon interactions with proton sources.") ;

  //This creates the singleton which deals with the scaling factor for IRB Mean free path for photo disintegration and pionproduction. To create the singleton, one needs the location of the xml file. Thus, GetInstance() is called here. All other calls can now use GetInstance() with the xml file name.
  //IRBzEvolutionFactory::GetInstance(string((aFileName)));
  //IRBzEvolutionFactory::GetInstance()->PlotAllModels();
}

TEveryThing::~TEveryThing() {

  if ( Basic ) delete Basic;
  if ( Univ  ) delete Univ;

}
