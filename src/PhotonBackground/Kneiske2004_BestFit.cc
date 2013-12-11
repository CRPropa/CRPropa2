#include <Kneiske2004_BestFit.h> 
#include <algorithm>

Kneiske2004_BestFit:: Kneiske2004_BestFit(std::string XMLFileName): TIRBzEvolutionModel(XMLFileName){
  std::string Path = DEFAULT_IRBzRedshiftEvol_Kneiske;
  _fMaxRedshift = 4;

  //Read Plots from paper page 8 fig. 4 into TGraphs (root). In this way they can be integrated.
  std::vector<TGraph*> TGraph_Vec;
  TGraph z_0_TG((Path+"PlotFrom0309141_z_0.csv").c_str(),"%lg,%lg");
  TGraph_Vec.push_back(&z_0_TG);
  TGraph z_0_2_TG((Path+"PlotFrom0309141_z_0_2.csv").c_str(),"%lg,%lg");
  TGraph_Vec.push_back(&z_0_2_TG);
  TGraph z_0_4_TG((Path+"PlotFrom0309141_z_0_4.csv").c_str(),"%lg,%lg");
  TGraph_Vec.push_back(&z_0_4_TG);
  TGraph z_0_6_TG((Path+"PlotFrom0309141_z_0_6.csv").c_str(),"%lg,%lg");
  TGraph_Vec.push_back(&z_0_6_TG);
  TGraph z_1_TG((Path+"PlotFrom0309141_z_1.csv").c_str(),"%lg,%lg");
  TGraph_Vec.push_back(&z_1_TG);
  TGraph z_2_TG((Path+"PlotFrom0309141_z_2.csv").c_str(),"%lg,%lg");
  TGraph_Vec.push_back(&z_2_TG);
  TGraph z_3_TG((Path+"PlotFrom0309141_z_3.csv").c_str(),"%lg,%lg");
  TGraph_Vec.push_back(&z_3_TG);
  TGraph z_4_TG((Path+"PlotFrom0309141_z_4.csv").c_str(),"%lg,%lg");
  TGraph_Vec.push_back(&z_4_TG);
  
  
  for(int i=0; i<TGraph_Vec.size(); i++){
    for(int j=0; j<TGraph_Vec[i]->GetN(); j++){
      double x,y;
      if(TGraph_Vec[i]->GetPoint(j, x, y)==-1) throw TCrpErr("Error in Kneiske2004_BestFit(): can't read point from TGraph."); 

      //convert wavelength from \mu m to m.
      x=x*1.e-6;

      //convert y axis from nW/(m^2*sr) to nW/(m^2*sr*m). Note, 1/m means per wavelength here!
      y=y/x;
      
      //convert nW=1.e-9*J/s=1.e-9*1.e19/1.602 eV/s=1.e-9*1.e19/1.602*1.e-9*GeV/s. Now, [y]=GeV/(m^2*sr*m*s).
      y=y*1.e-9/(1.602e-19)*1.e-9;

      //convert the 1/m (which means per wavelength) factor to 1/GeV (which means per photon energy). 
      //Thus, multiply with x/(h*\nu)=x*x/(h*c) => [y]=GeV/(m^2*sr*GeV*s)
      y=y*(x*x/(4.136e-24*3.e8));

      //Now calculate the energy per unit volume (s. the argumentation as in Carol and Ostlie (An introdutcion to modern astrophysics, Addison Wesley, 1996)).
      y = y * 4 * M_PI / 3.e8;

      //convert wavelength x from m to photon energy in GeV. 
      x=4.136e-24/(x/3.e8);
      
      //Divide by photon energy to yield number of photons #/(m^3*GeV)
      y=y/x;

      //In units h_bar=c=1 the units are GeV^2. Thus, multiply with (h_bar*c)^3 which has units of (GeV*m)^3.
      y=y * pow( 6.528e-25 * 3.e8 , 3 );
      
      TGraph_Vec[i]->SetPoint( j, /*log10*/ (x) , y );

    }
  }
  


  TCanvas PowerSpectraCan("PowerSpectraCan","PowerSpectraCan",1);
  z_0_TG.Draw("AC");
  z_0_TG.SetLineStyle(0);
  z_0_TG.GetYaxis()->SetTitle("photon number density n(#epsilon) / GeV^2");
  z_0_TG.GetXaxis()->SetTitle("photon energy #epsilon / GeV");
    
  z_0_2_TG.Draw("C");
  z_0_2_TG.SetLineStyle(1);

  z_0_4_TG.Draw("C");
  z_0_4_TG.SetLineStyle(2);  

  z_0_6_TG.Draw("C");
  z_0_6_TG.SetLineStyle(3);

  z_1_TG.Draw("C");
  z_1_TG.SetLineStyle(4);
    
  z_2_TG.Draw("C");
  z_2_TG.SetLineStyle(5);
  
  z_3_TG.Draw("C");
  z_3_TG.SetLineStyle(6);

  z_4_TG.Draw("C");
  z_4_TG.SetLineStyle(7);
  PowerSpectraCan.SetLogx();
  PowerSpectraCan.SetLogy();
  PowerSpectraCan.SaveAs("PowerSpectraCan_Kneiske2004_BestFit.root");


  TCanvas PowerSpectraCan_fineBin("PowerSpectraCan_fineBin","PowerSpectraCan_fineBin",1);
  int NIntPoints=10000;
  //copy TGraph with more points to find the intersection point of IRB and CMB numerically.
  std::vector<TGraph*> TGraph_Vec_fineBin;
  std::vector<double> CMB_X_Vec, CMB_Y_Vec;
  std::vector<double> Maxima_Vec;
  CMB CMB;
   for(int i=0; i<TGraph_Vec.size(); i++){
     double EMin, EMax, bufferY;
     TGraph_Vec[i]->GetPoint(TGraph_Vec[i]->GetN()-1, EMin, bufferY);
     double logEMin = log10(EMin);
     TGraph_Vec[i]->GetPoint(0, EMax, bufferY);
     double logEMax = log10(EMax);

     //use log binning
     std::vector<double> X_Vec, Y_Vec;
     for(double ActLogE=logEMin; 
	 ActLogE<logEMax; 
	 ActLogE+=fabs(logEMax-logEMin)/((double) NIntPoints)){

       if (i==1) {
	 CMB_Y_Vec.push_back(CMB.GetPhotonDensity(0.,0.,0.,0.,pow(10,ActLogE))); 
	 CMB_X_Vec.push_back(pow(10,ActLogE)); 
       }

       double ActYVal= TGraph_Vec[i]->Eval(pow(10.,ActLogE),0,"");
       //only use point for which the IRB dominates the CMB photons
       if(CMB.GetPhotonDensity(0.,0.,0.,0.,pow(10,ActLogE))>ActYVal) continue;

       X_Vec.push_back(pow(10.,ActLogE));
       Y_Vec.push_back(ActYVal);
     }
     TGraph* ActTG= new TGraph(X_Vec.size(), &(X_Vec[0]), &(Y_Vec[0]));
     TGraph_Vec_fineBin.push_back(ActTG);
     //find maximum
     Maxima_Vec.push_back(*std::max_element(Y_Vec.begin(), Y_Vec.end()));
     if(i==0) ActTG->Draw("A*l");
     else ActTG->Draw("*l");
   }
   TGraph* CMBTG= new TGraph(CMB_X_Vec.size(), &(CMB_X_Vec[0]), &(CMB_Y_Vec[0]));
   CMBTG->Draw("*l");
   //PowerSpectraCan_fineBin.SaveAs("PowerSpectraCan_Kneiske2004_BestFit_fineBin.root");
    
   
  std::vector<double> x, y;
  
  //normalization depends now on redshift. Old version with redshift-independent normalization is commented out


//  double norm = TGraph_Vec_fineBin[0]->Integral(0,TGraph_Vec_fineBin[0]->GetN())*pow(1+0.,3.);

  double norm0 = TGraph_Vec_fineBin[0]->Integral(0,TGraph_Vec_fineBin[0]->GetN())*pow(1+0.,3.);
  double norm02 = TGraph_Vec_fineBin[0]->Integral(0,TGraph_Vec_fineBin[1]->GetN())*pow(1+0.,3.);
  double norm04 = TGraph_Vec_fineBin[0]->Integral(0,TGraph_Vec_fineBin[2]->GetN())*pow(1+0.,3.);
  double norm06 = TGraph_Vec_fineBin[0]->Integral(0,TGraph_Vec_fineBin[3]->GetN())*pow(1+0.,3.);
  double norm1 = TGraph_Vec_fineBin[0]->Integral(0,TGraph_Vec_fineBin[4]->GetN())*pow(1+0.,3.);
  double norm2 = TGraph_Vec_fineBin[0]->Integral(0,TGraph_Vec_fineBin[5]->GetN())*pow(1+0.,3.);
  double norm3 = TGraph_Vec_fineBin[0]->Integral(0,TGraph_Vec_fineBin[6]->GetN())*pow(1+0.,3.);
  double norm4 = TGraph_Vec_fineBin[0]->Integral(0,TGraph_Vec_fineBin[7]->GetN())*pow(1+0.,3.);

//  x.push_back(0.); y.push_back( TGraph_Vec_fineBin[0]->Integral(0,TGraph_Vec_fineBin[0]->GetN())*pow(1+0.,3.)/norm);
//  x.push_back(0.2); y.push_back( TGraph_Vec_fineBin[1]->Integral(0,TGraph_Vec_fineBin[1]->GetN())*pow(1+0.2,3.)/norm);
//  x.push_back(0.4); y.push_back( TGraph_Vec_fineBin[2]->Integral(0,TGraph_Vec_fineBin[2]->GetN())*pow(1+0.4,3.)/norm);
//  x.push_back(0.6); y.push_back( TGraph_Vec_fineBin[3]->Integral(0,TGraph_Vec_fineBin[3]->GetN())*pow(1+0.6,3.)/norm);
//  x.push_back(1); y.push_back( TGraph_Vec_fineBin[4]->Integral(0,TGraph_Vec_fineBin[4]->GetN())*pow(1+1.,3.)/norm);
//  x.push_back(2); y.push_back( TGraph_Vec_fineBin[5]->Integral(0,TGraph_Vec_fineBin[5]->GetN())*pow(1+2.,3.)/norm);
//  x.push_back(3); y.push_back( TGraph_Vec_fineBin[6]->Integral(0,TGraph_Vec_fineBin[6]->GetN())*pow(1+3.,3.)/norm);
//  x.push_back(4); y.push_back( TGraph_Vec_fineBin[7]->Integral(0,TGraph_Vec_fineBin[7]->GetN())*pow(1+4.,3.)/norm);

  x.push_back(0.); y.push_back( TGraph_Vec_fineBin[0]->Integral(0,TGraph_Vec_fineBin[0]->GetN())*pow(1+0.,3.)/norm0);
  x.push_back(0.2); y.push_back( TGraph_Vec_fineBin[1]->Integral(0,TGraph_Vec_fineBin[1]->GetN())*pow(1+0.2,3.)/norm02);
  x.push_back(0.4); y.push_back( TGraph_Vec_fineBin[2]->Integral(0,TGraph_Vec_fineBin[2]->GetN())*pow(1+0.4,3.)/norm04);
  x.push_back(0.6); y.push_back( TGraph_Vec_fineBin[3]->Integral(0,TGraph_Vec_fineBin[3]->GetN())*pow(1+0.6,3.)/norm06);
  x.push_back(1); y.push_back( TGraph_Vec_fineBin[4]->Integral(0,TGraph_Vec_fineBin[4]->GetN())*pow(1+1.,3.)/norm1);
  x.push_back(2); y.push_back( TGraph_Vec_fineBin[5]->Integral(0,TGraph_Vec_fineBin[5]->GetN())*pow(1+2.,3.)/norm2);
  x.push_back(3); y.push_back( TGraph_Vec_fineBin[6]->Integral(0,TGraph_Vec_fineBin[6]->GetN())*pow(1+3.,3.)/norm3);
  x.push_back(4); y.push_back( TGraph_Vec_fineBin[7]->Integral(0,TGraph_Vec_fineBin[7]->GetN())*pow(1+4.,3.)/norm4);
  
  /*
  double norm = Maxima_Vec[0]*pow(1+0.,3.);
  x.push_back(0.); y.push_back( Maxima_Vec[0]*pow(1+0.,3.)/norm); 
  x.push_back(0.2); y.push_back( Maxima_Vec[1]*pow(1+0.2,3.)/norm); 
  x.push_back(0.4); y.push_back( Maxima_Vec[2]*pow(1+0.4,3.)/norm); 
  x.push_back(0.6); y.push_back( Maxima_Vec[3]*pow(1+0.6,3.)/norm); 
  x.push_back(1); y.push_back( Maxima_Vec[4]*pow(1+1.,3.)/norm); 
  x.push_back(2); y.push_back( Maxima_Vec[5]*pow(1+2.,3.)/norm); 
  x.push_back(3); y.push_back( Maxima_Vec[6]*pow(1+3.,3.)/norm); 
  x.push_back(4); y.push_back( Maxima_Vec[7]*pow(1+4.,3.)/norm); 
  */
  TCanvas IRBScaleFactorCan("IRBScaleFactorCan","IRBScaleFactorCan",1);
  fIRBScaleFactorTG = new TGraph(x.size(),&(x[0]),&(y[0]));
  fIRBScaleFactorTG->Draw("ACp");
  IRBScaleFactorCan.SaveAs("IRBScaleFactorCan.root");
}

double Kneiske2004_BestFit::GetScalingFactor(double redshift){
  if(redshift>_fMaxRedshift) return 0.;  

  return fIRBScaleFactorTG->Eval(redshift,0,"");
}
