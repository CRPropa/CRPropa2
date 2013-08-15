/**
@file    CRPropa.cc
 @author  Eric Armengaud, armengau@in2p3.fr
 @brief   Main file for UHECR propagation code
*/

#include "CRPropa.h" 
#include <exception>

double interactt, deflectt, photodt, pionprodt,  pairpt, decayt;
long pprod, bbarprod, nprod;

int main (int argc, char **argv) {

  try{
    if ( argc < 2 ) throw TCrpErr("No XML file specified as arg.") ;

    // Build basic parameters and environment    
    TEveryThing lAll(argv[1]);
    TBasicParam *lpBasic = lAll.Basic;
    TUniverse   *lpUniv  = lAll.Univ;


    pprod=0;
    bbarprod=0;
    nprod=0;

    //TODO: check if this can be done in a different way! 
    if(
       lpUniv->Observers()->Type() == OBSERVER_SMALLSPHERE 
       && 
       (lpUniv->IntegratorMinTimeStep() * c_light / Mpc >  lpUniv->Observers()->Radius()*OBS_SECURITY_RATIO/Mpc )
       ) 
      TCrpErr("In small oberserver mode MinStep_Mpc should be smaller than Radius_Mpc*OBS_SECURITY_RATIO.");
 
    QUEUE<TParticle*> lParts ;          // List of particles to propagate
    TList1DPhotons lPhotons1D(lpUniv) ; // List of photons IN THE 1D CASE
    TList1DPhotons lPhotons1DPrimary(lpUniv, "PRIMARY") ; // List of photons IN THE 1D CASE
    unsigned long lN = 0;
    std::cout<<"Start nuclei propagation Progress : ";
    while ( lN < lpBasic->N() ) {
      if (lpUniv->Sources()->IsPhotonSources()) { // Photon injection
	lPhotons1DPrimary.AddPhoton(lpUniv) ;
	lN++ ;
      } else {
	lParts.clear(); // just to be sure.
	lParts.push_back( new TNucleus(lpUniv, &lPhotons1D ));

	if(fmod((double) lN*100./lpBasic->N(), 10.) <= std::numeric_limits<double>::min()){  std::cout<<"# "<<(int) (lN*100./lpBasic->N())<<"% #"<<std::endl;}

	lN++;
	while ( lParts.size() ) {
	  QUEUE<TParticle*>* lNewSet;
	  //Propagation
	  lNewSet = lParts.back()->Propagate(lpBasic) ;
	  //Delete propagated particle
	  delete lParts.back() ;
	  lParts.pop_back();
	  //Adding to main set all new particles
	  lParts.insert(lParts.end(),lNewSet->begin(),lNewSet->end());
	  lNewSet->clear();
	}
      }
    }
    
    std::cout<<"# 100% \n\t-> trajectories done!"<<std::endl;

    std::cout<<"Start propagation of secondaries if selected . . . this may take some time."<<std::endl;

    lPhotons1DPrimary.Propagate(lpUniv,lpBasic) ; // Propagate primary photons
    lPhotons1D.Propagate(lpUniv,lpBasic) ; // IN THE 1D CASE

    cerr << argv[0] << " : Normal end." << endl;
   
#ifdef DEBUG_OUTPUT
    std::cout << "interact : " << interactt / CLOCKS_PER_SEC << std::endl;
    std::cout << "deflect : " << deflectt / CLOCKS_PER_SEC << std::endl;
    std::cout << "photod : " << photodt / CLOCKS_PER_SEC << std::endl;
    std::cout << "pionprod : " << pionprodt / CLOCKS_PER_SEC << std::endl;
    std::cout << "pairp : " << pairpt / CLOCKS_PER_SEC << std::endl;
    std::cout << "decayt : " << decayt / CLOCKS_PER_SEC << std::endl;
    std::cout << "pprod : " << pprod << std::endl;
    std::cout << "bbarprod : " << bbarprod << std::endl;
    std::cout << "nprod : " << nprod << std::endl;

#endif
    
#ifdef NILS 
    vector< unsigned long int > NucleiAvailable 
      = TabulatedTALYSCrossSection::GetInstance()->GetListOfInputNuclei();
    vector< double > MassVec;
    vector< double > ChargeVec;
    for(int i=0; i<NucleiAvailable.size(); i++){
      std::cout<<"NucleiAvailable[]="<<NucleiAvailable[i]<<std::endl;
      vector< unsigned long int > MassAndCharg
	=TabulatedTALYSCrossSection::GetInstance()->GetMassAndChargeNuclId(NucleiAvailable[i]);
      ChargeVec.push_back(MassAndCharg[1]+0.5);
      MassVec.push_back(MassAndCharg[0]+0.5);      
      std::cout<<"M="<<MassAndCharg[0]<<"\tZ="<<MassAndCharg[1]<<std::endl;
    }
    MassChargePathCan->cd();
    MassChargePathTH2F->Draw("colz");    
    TGraph* MassChargePathAddOnTG = new TGraph(ChargeVec.size(),
					       &(ChargeVec[0]),
					       &(MassVec[0]));
    MassChargePathAddOnTG->Draw("*");
    MassChargePathAddOnTG->SetFillColor(0);
    MassChargePathAddOnTG->SetMarkerStyle(4);	

    //Rachen
    vector< double > MassVec_Rachen;
    vector< double > ChargeVec_Rachen;
    MassVec_Rachen.push_back(2+0.5);ChargeVec_Rachen.push_back(1+0.5);
    MassVec_Rachen.push_back(3+0.5);ChargeVec_Rachen.push_back(1+0.5);
    MassVec_Rachen.push_back(3+0.5);ChargeVec_Rachen.push_back(2+0.5);
    MassVec_Rachen.push_back(4+0.5);ChargeVec_Rachen.push_back(2+0.5);
    TGraph* MassChargePathRachenTG = new TGraph(ChargeVec_Rachen.size(),
					       &(ChargeVec_Rachen[0]),
					       &(MassVec_Rachen[0]));
    //MassChargePathTH2F->Draw("colzsame");    
    MassChargePathRachenTG->Draw("*");
    MassChargePathRachenTG->SetFillColor(0);
    MassChargePathRachenTG->SetMarkerStyle(20);
    
    //Our add ons
    vector< double > MassVec_OurAdds;
    vector< double > ChargeVec_OurAdds;
    MassVec_OurAdds.push_back(7+0.5);ChargeVec_OurAdds.push_back(3+0.5);
    MassVec_OurAdds.push_back(9+0.5);ChargeVec_OurAdds.push_back(4+0.5);
    TGraph* MassChargePathOurAddsTG = new TGraph(ChargeVec_OurAdds.size(),
						 &(ChargeVec_OurAdds[0]),
						 &(MassVec_OurAdds[0]));
    MassChargePathOurAddsTG->Draw("*");
    MassChargePathOurAddsTG->SetMarkerStyle(21);
    MassChargePathOurAddsTG->SetFillColor(0);
    MassChargePathOurAddsTG->SetMarkerSize(1.2);
    
    //Geant 4 parametrization
    vector< double > Geant4MassVec_OurAdds;
    vector< double > Geant4ChargeVec_OurAdds;

    //lithium
    //Geant4MassVec_OurAdds.push_back(7+0.5);  Geant4ChargeVec_OurAdds.push_back(3+0.5); 
    Geant4MassVec_OurAdds.push_back(8+0.5);  Geant4ChargeVec_OurAdds.push_back(3+0.5); 
    Geant4MassVec_OurAdds.push_back(9+0.5);  Geant4ChargeVec_OurAdds.push_back(3+0.5); 
 
    //berylium
    Geant4MassVec_OurAdds.push_back(7+0.5);  Geant4ChargeVec_OurAdds.push_back(4+0.5); 
    //Geant4MassVec_OurAdds.push_back(9+0.5);  Geant4ChargeVec_OurAdds.push_back(4+0.5); 
    Geant4MassVec_OurAdds.push_back(10+0.5);  Geant4ChargeVec_OurAdds.push_back(4+0.5); 
    Geant4MassVec_OurAdds.push_back(11+0.5);  Geant4ChargeVec_OurAdds.push_back(4+0.5); 
    
    //bor
    Geant4MassVec_OurAdds.push_back(8+0.5);  Geant4ChargeVec_OurAdds.push_back(5+0.5); 
    Geant4MassVec_OurAdds.push_back(10+0.5);  Geant4ChargeVec_OurAdds.push_back(5+0.5); 
    Geant4MassVec_OurAdds.push_back(11+0.5);  Geant4ChargeVec_OurAdds.push_back(5+0.5); 
    
    //carbon
    Geant4MassVec_OurAdds.push_back(9+0.5);  Geant4ChargeVec_OurAdds.push_back(6+0.5); 
    Geant4MassVec_OurAdds.push_back(10+0.5);  Geant4ChargeVec_OurAdds.push_back(6+0.5); 
    Geant4MassVec_OurAdds.push_back(11+0.5);  Geant4ChargeVec_OurAdds.push_back(6+0.5); 
   
    TGraph* Geant4MassChargePathOurAddsTG = new TGraph(Geant4ChargeVec_OurAdds.size(),
						       &(Geant4ChargeVec_OurAdds[0]),
						       &(Geant4MassVec_OurAdds[0]));
    Geant4MassChargePathOurAddsTG->Draw("*");
    Geant4MassChargePathOurAddsTG->SetMarkerStyle(25);
    Geant4MassChargePathOurAddsTG->SetFillColor(0);
    Geant4MassChargePathOurAddsTG->SetMarkerSize(1.6);

    TLine* TALYSRel   = new TLine(0,12,9,12);
    TALYSRel->SetLineColor(8);
    TALYSRel->SetLineWidth(2);
    TALYSRel->Draw();
    TLine* TALYSUsage = new TLine(0,9,7,9);
    TALYSUsage->SetLineColor(5);
    TALYSUsage->SetLineWidth(2);
    TALYSUsage->Draw();
    TLine* TALYSNoGo = new TLine(0,6,4,6);
    TALYSNoGo->SetLineColor(2);
    TALYSNoGo->SetLineWidth(2);
    TALYSNoGo->Draw();
    TF1* NLarger2Line = new TF1("NL2","2+(x)+1",0.,27.);
    //NLarger2Line->SetLineColor(2);
    //NLarger2Line->SetLineStyle(2);
    //NLarger2Line->Draw("same"); 
    std::vector<double> vectorX;
    std::vector<double> vectorY;
    for(double d=4.; d<27; d+=0.001){
      vectorX.push_back(d);
      vectorY.push_back((int) NLarger2Line->Eval(d));
    }
    
    for(int n=4.; n<=27; n++){
      double lower = (int) NLarger2Line->Eval(n-0.5);
      double upper = (int) NLarger2Line->Eval(n+0.5);
      for(double y=lower; y<upper; y+=0.001){
	vectorX.push_back(n);
	vectorY.push_back(y);
      }
    }

    TGraph*NLarger2LineTG = new TGraph(vectorX.size(), &(vectorX[0]), &(vectorY[0]));
    NLarger2LineTG->SetMarkerColor(2);
    NLarger2LineTG->SetMarkerStyle(20);
    NLarger2LineTG->SetMarkerSize(0.3);
    NLarger2LineTG->Draw("p");

    TLegend *legend = new TLegend(.75,.80,.95,.95);
    legend->AddEntry(TALYSRel,"TALYS realiablity Kroning et al.");
    legend->AddEntry(TALYSUsage,"TALYS \"reliability\" Allard et al.");
    legend->AddEntry(TALYSNoGo,"TALYS abortion line 1 (A>5 or N>2)");
    //legend->AddEntry(NLarger2LineTG,"TALYS abortion line 2 (N>2)");
    legend->AddEntry(MassChargePathAddOnTG,"TALYS cross section available");
    legend->AddEntry(MassChargePathRachenTG,"cross sections Rachenet al.");
    legend->AddEntry(Geant4MassChargePathOurAddsTG,"cross sections as in Geant4");
    legend->AddEntry(MassChargePathOurAddsTG,"available measured cross section");
    legend->Draw();
    MassChargePathCan->SaveAs("MassChargePathCan.root");

    MassChargeVSChannelCan->cd(); 
    MassChargeVSChannelTH2F->Draw("colz"); 
    MassChargeVSChannelCan->SaveAs("MassChargeVSChannelCan.root");
    MFPvsEnCan->cd();
    MFPvsEnTH2F->Draw("colz");
    MFPvsEnTH2F->ProfileX()->Draw("same");
    TInteractionData *lpInt = lAll.Univ->InteractionData() ;
    std::vector<double> XVal, YVal, YVal_PiP, YVal_PiPAndPD, YVal_PD_Excl;
    for(double i=MFPvsEnTH2F->GetXaxis()->GetBinCenter(0); 
	i<=MFPvsEnTH2F->GetXaxis()->GetBinCenter(MFPvsEnTH2F->GetNbinsX()); 
	i+=(MFPvsEnTH2F->GetXaxis()->GetBinCenter(MFPvsEnTH2F->GetNbinsX())-MFPvsEnTH2F->GetXaxis()->GetBinCenter(0))/50.){
      double log10_gamma = i;//MFPvsEnTH2F->GetXaxis()->GetBinCenter(i);
      XVal.push_back(log10_gamma);

      double ActMFP = lpInt->GetPDTimeStep(26*1000 + 56,
			 0,
			 pow(10,log10_gamma)*56,     //GeV
			 56,             //GeV
			 0.);
      if(ActMFP!=0){ ActMFP = 1/(ActMFP*(Mpc * inv_c_light));}
      else{ActMFP=1.e99;}
      //Convert MFP to expectation value E[log10(MFP)] in the logarithmic binning
      std::ostringstream Formula;
      Formula<<"log10(x)"<<"*exp(-x/"<<ActMFP<<")/"<<ActMFP;
      std::cout<<"string: "<<Formula.str().c_str()<<std::endl;
      TF1* convertMFPToLog10Binning = new TF1("convertMFPToLog10Binning",Formula.str().c_str(),0, 20*ActMFP);
      double ActMFP_logBinning = convertMFPToLog10Binning->Integral(0, 20*ActMFP);
      //YVal.push_back(log10(ActMFP));
      YVal.push_back(ActMFP_logBinning);

      TNucleus* nucleus = new TNucleus(lAll.Univ, 
				   0., 
				   pow(10,log10_gamma)*56 * 1000 * inv_c_light, 
				   0., 
				   pow(10,log10_gamma)*56 * 1000 * inv_c_light, 
				   3400., 
				   56,
				   26,
				   &lPhotons1DPrimary, 
				   26056
				   );

      double ActMFP_PiP = (   1/nucleus->PionProtonMFP( 1) + 1/nucleus->PionNeutronMFP( 1)
			      +
			      1/nucleus->PionProtonMFP( 2) + 1/nucleus->PionNeutronMFP( 2)  )    *    (Mpc * inv_c_light);
      ActMFP_PiP=1/(ActMFP_PiP);
      
      std::ostringstream Formula_PiP;
      Formula_PiP<<"log10(x)"<<"*exp(-x/"<<ActMFP_PiP<<")/"<<ActMFP_PiP;
      std::cout<<"string: "<<Formula_PiP.str().c_str()<<std::endl;
      TF1* convertMFPToLog10Binning_PiP = new TF1("convertMFPToLog10Binning_PiP",Formula_PiP.str().c_str(),0, 20*ActMFP_PiP);
      double ActMFP_PiP_logBinning = convertMFPToLog10Binning_PiP->Integral(0, 20*ActMFP_PiP);
      //YVal.push_back(log10(ActMFP));
      YVal_PiP.push_back(ActMFP_PiP_logBinning);

      double ActMFP_PiPAndPD=1. / ( 1./ActMFP + 1./ActMFP_PiP ) ;
      std::ostringstream Formula_PiPAndPD;
      Formula_PiPAndPD<<"log10(x)"<<"*exp(-x/"<<ActMFP_PiPAndPD<<")/"<<ActMFP_PiPAndPD;
      std::cout<<"string: "<<Formula_PiPAndPD.str().c_str()<<std::endl;
      TF1* convertMFPToLog10Binning_PiPAndPD = new TF1("convertMFPToLog10Binning_PiPAndPD",Formula_PiPAndPD.str().c_str(),0, 20*ActMFP_PiPAndPD);
      double ActMFP_PiPAndPD_logBinning = convertMFPToLog10Binning_PiPAndPD->Integral(0, 20*ActMFP_PiPAndPD);

      YVal_PiPAndPD.push_back(ActMFP_PiPAndPD_logBinning);
      
      double MFP_PD_Excl_CMB=TabulatedTALYSMeanFreePath_CMB::GetInstance()->GetExclusiveChannelY(26056,
												 10000,
												 log10_gamma);
      ostringstream buffer;
      buffer<<DEFAULT_TALYSTabMeanFreePath_DIR;
      TalysDirectory = buffer.str();
      double MFP_PD_Excl_IRB=TabulatedTALYSMeanFreePath_IRB::GetInstance()->GetExclusiveChannelY(26056,
												 10000,
												 log10_gamma);
      YVal_PD_Excl.push_back( log10( 1/(MFP_PD_Excl_CMB+MFP_PD_Excl_IRB) ) );

      std::cout<<"log10_gamma="<<log10_gamma
	       <<"\t En="<<pow(10,log10_gamma)*56*1000
	       <<"\t MFP-PD  ="<<ActMFP
	       <<"\t MFP-PiP ="<<ActMFP_PiP
	       <<"\t MFP excl="<<1/(MFP_PD_Excl_CMB+MFP_PD_Excl_IRB)
	       <<std::endl;
    }
    TGraph* Overlay     = new TGraph(XVal.size(), &(XVal[0]),  &(YVal[0]));
    TGraph* Overlay_PiP = new TGraph(XVal.size(), &(XVal[0]),  &(YVal_PiP[0]));
    TGraph* Overlay_PiPAndPD = new TGraph(XVal.size(), &(XVal[0]),  &(YVal_PiPAndPD[0]));
    MFPvsEnCan->cd();
    Overlay->Draw("l");
    Overlay_PiP->Draw("l");
    Overlay_PiPAndPD->Draw("l");
    MFPvsEnCan->SaveAs("MFPvsEnCan.root");

    
    TGraph* Overlay_PDExcl = new TGraph(XVal.size(), &(XVal[0]),  &(YVal_PD_Excl[0]));
    MFPvsEn_ExclCan->cd();
    MFPvsEn_ExclTH2F->Draw("colz");
    MFPvsEn_ExclTH2F->ProfileX()->Draw("same");
    Overlay_PDExcl->Draw("l");
    MFPvsEn_ExclCan->SaveAs("MFPvsEn_ExclCan.root");

    //  //Small sanity check
    //     NucleiAvailable .clear();
    //     NucleiAvailable 
    //       = TabulatedTALYSCrossSection::GetInstance()->GetListOfInputNuclei();
    //     MassVec.clear();
    //     ChargeVec.clear();
    //     for(int i=0; i<NucleiAvailable.size(); i++){
    //       std::cout<<"\n*****NucleiAvailable[]="<<NucleiAvailable[i]<<std::endl;
    //       vector< unsigned long int > MassAndCharg
    // 	=TabulatedTALYSCrossSection::GetInstance()->GetMassAndChargeNuclId(NucleiAvailable[i]);
    //       ChargeVec.push_back(MassAndCharg[1]+0.5);
    //       MassVec.push_back(MassAndCharg[0]+0.5);      
    //       std::cout<<"M="<<MassAndCharg[0]<<"\tZ="<<MassAndCharg[1]<<std::endl;
    //       vector< unsigned long int > XSAvailable 
    // 	= TabulatedTALYSCrossSection::GetInstance()->GetListOfAllExclusiveChannels(NucleiAvailable[i]);
    
    //       for(int m=0; m<XSAvailable.size(); m++){
    // 	 std::vector< int > DMassAndDChargeExcl 
    // 	   = TabulatedTALYSCrossSection::GetInstance()->GetDMassAndDChargeExclId (XSAvailable[m]);
    // 	 std::cout<<"xs id="<<XSAvailable[m]
    // 		  <<"\t DM="<< DMassAndDChargeExcl[0]
    // 		  <<"\t DZ="<< DMassAndDChargeExcl[1]
    // 	   <<std::endl;
    // 	 if(MassAndCharg[0]<DMassAndDChargeExcl[0]){
    // 	   std::cout<<"AHHHHH nucleusId="<<NucleiAvailable[i]<<" xsid="<<XSAvailable[m]<<std::endl;
    // 	   //exit(-1);
    // 	 }	 
    //       }
    
    //     }
    
#endif


    return 0;
    
  } catch (exception& e ) {
      cout << e.what() << endl;
    cerr << "Caught exception but, unfortunately, do not know what to do with."
	 << " Exiting." << endl;
    exit(ERR_MAIN);
  }
  
}
