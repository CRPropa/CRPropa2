/**
   @file    sophiainteractions.cc
   @author  Eric Armengaud, armengau@in2p3.fr
   @brief   Implementation of the TSophiaInteractions class. See the .h file
*/

#include "sophiainteractions.h"
#include "nucleusdb.h"
#include "TabulatedTALYSAveragedCrossSection.h"
#include <TALYSMeanFreePathAvrgd.h>
#include <TabulatedTALYSMeanFreePath_CMB.h>
#include <TabulatedTALYSMeanFreePath_IRB.h>
#include <iostream>
#include <PhotonBackground.h>
#include <CMB.h>
#include <IR.h>
#include <variableinfrared.h>
#include <CRPropa.h>
#include <algorithm>

std::string TalysDirectory;

//#include "nucleus.cc"

//This is needed for the random_shuffle() funtion from algorithm.h (C++ stdl). 
// random generator function:
ptrdiff_t myrandom (ptrdiff_t i) { return RandFlat::shootInt(i);}
// pointer object to it:
ptrdiff_t (*p_myrandom)(ptrdiff_t) = myrandom;

TSophiaInteractions::TSophiaInteractions(const char* aFileName, TIRBzEvolutionModel* _lpModel) : TXmlParam(aFileName)
{
  SetType(INTERACTION_SOPHIA);

    _fpEvolModel = _lpModel;
    
  string lIntDir ;
  TiXmlElement* lpXmlInter = XmlExtract().GetElement("Interactions") ;
  
  // Interaction directory
  TiXmlElement* lpXmlDir = lpXmlInter->FirstChildElement("Directory") ;
  if ( ! lpXmlDir ) {
    //throw TXmlErr("Incorrect interaction directory");
    // if undef then use default one
    lIntDir = DEFAULT_SOPHIA_DIR ;
  }
  else {
    TiXmlNode* lpXmlDirName = lpXmlDir->FirstChild() ;
    if ( (! lpXmlDirName) || (! lpXmlDirName->ToText()) )
      throw TXmlErr("Incorrect interaction directory name");
    lIntDir = lpXmlDirName->Value() ;
  }
  
  // Maximum time step for interactions
  TiXmlElement* lpXmlStep = lpXmlInter->FirstChildElement("MaxStep_Mpc") ;
  if ( ! lpXmlStep ) 
    throw TXmlErr("No MaxStep_Mpc for interactions") ;
  const char* lStep = lpXmlStep->Attribute("value",&_fInteractionTimeStep) ;
  if ( ! lStep) 
    throw TXmlErr("Error getting value of interaction step in xml file");
  _fInteractionTimeStep *= Mpc * inv_c_light;
  
  // Interactions and secondary propagation on or off
  TiXmlElement* lpXmlCheckRedshift = lpXmlInter->FirstChildElement("NoRedshift") ;
  if (!lpXmlCheckRedshift) {
    _fRedshiftFlag = 1 ;
  } else _fRedshiftFlag = 0 ;
  TiXmlElement* lpXmlCheckNoRedshiftEnergyLoss = lpXmlInter->FirstChildElement("NoRedshiftEnergyLoss") ;
  if (!lpXmlCheckNoRedshiftEnergyLoss) {
    _fNoRedshiftEnergyLossFlag = 1 ;
  } else _fNoRedshiftEnergyLossFlag = 0 ;
  TiXmlElement* lpXmlCheckPhoton = lpXmlInter->FirstChildElement("SecondaryPhotons") ;
  if ( lpXmlCheckPhoton ) {
    _fSecPhotonFlag = 1 ;
  } else _fSecPhotonFlag = 0 ;
  TiXmlElement* lpXmlCheckPairProdPhoton = lpXmlInter->FirstChildElement("SecondaryPairProdPhotons") ;
  _fKelnerPairProdFlag=1;
  if ( lpXmlCheckPairProdPhoton ) {
    _fSecPairProdPhotonFlag = 1 ;
    const char* lProba = lpXmlCheckPairProdPhoton->Attribute("proba",&_fSecPairProdPhotonProba) ;
    if ( ! lProba ) _fSecPairProdPhotonProba = 1 ;
    if ( _fSecPairProdPhotonProba < 0 || _fSecPairProdPhotonProba > 1 )
      throw TXmlErr("Probability not in [0;1] for SecondaryPairProdPhotons");
    const char* toto = lpXmlCheckPairProdPhoton->Attribute("type");
    if (toto) {
      string lKelnerFlag = lpXmlCheckPairProdPhoton->Attribute("type") ;
      if (lKelnerFlag=="Power Law") _fKelnerPairProdFlag=0;
    }
  } else _fSecPairProdPhotonFlag = 0 ;
  TiXmlElement* lpXmlPPEps = lpXmlInter->FirstChildElement("PairProd_Eps") ;
  if (lpXmlPPEps ) {
    const char* lPPEps = lpXmlPPEps->Attribute("value",&_fPairProdEps) ;

    if ( !lPPEps ) _fPairProdEps = .01 ;
    if ( _fPairProdEps < 0 || _fPairProdEps > .01 )
      std::cerr << "Warning PairProd_Eps not in [0;.01]" << std::cerr;
  } else _fPairProdEps  = 0.01 ;
  TiXmlElement* lpXmlCheckNu = lpXmlInter->FirstChildElement("SecondaryNeutrinos") ;
  if ( lpXmlCheckNu ) {
    _fSecNuFlag = 1 ;
  } else _fSecNuFlag = 0 ;
  TiXmlElement* lpXmlCheckPairP = lpXmlInter->FirstChildElement("NoPairProd") ;
  if ( !lpXmlCheckPairP ) {
    _fPairProdFlag = 1 ;
  } else _fPairProdFlag = 0 ;
  TiXmlElement* lpXmlCheckPionP = lpXmlInter->FirstChildElement("NoPionProd") ;
  if ( !lpXmlCheckPionP ) {
    _fPionProdFlag = 1 ;
  } else _fPionProdFlag = 0 ;
  
  TiXmlElement* lpXmlCheckIRPionP = lpXmlInter->FirstChildElement("NoIRO") ;
  if ( !lpXmlCheckIRPionP ) {

    lpXmlCheckIRPionP = lpXmlInter->FirstChildElement("NoIRPionProd") ;
    if( !lpXmlCheckIRPionP) _fIRPionProdFlag=1; else _fIRPionProdFlag=0;
  
    TiXmlElement* lpXmlCheckIRPhotodisintegrationFlag = lpXmlInter->FirstChildElement("NoIRPhotodisintegration") ;
    if ( !lpXmlCheckIRPhotodisintegrationFlag )  _fIRPhotodisintegrationFlag = 1 ;  else _fIRPhotodisintegrationFlag = 0 ;
 
  } else _fIRPionProdFlag = 0 ;

  TiXmlElement* lpXmlCheckCMBPhotodisintegrationFlag = lpXmlInter->FirstChildElement("NoCMBPhotodisintegration") ;
  if ( !lpXmlCheckCMBPhotodisintegrationFlag ) {
    _fCMBPhotodisintegrationFlag = 1 ;
  } else _fCMBPhotodisintegrationFlag = 0 ;

  TiXmlElement* lpXmlCheckPhotodisintegration = lpXmlInter->FirstChildElement("NoPhotodisintegration") ;
  if ( !lpXmlCheckPhotodisintegration ) {
    _fPhotodisintegrationFlag = 1 ;
  } else _fPhotodisintegrationFlag = 0 ;

  TVector3D Positon; Positon.setX(0.); Positon.setY(0.); Positon.setZ(0.);
  TiXmlElement* lpXmlCheckPhotodisintegration_VarIRB = lpXmlInter->FirstChildElement("Photodisintegration_VarIRB") ;
  if ( lpXmlCheckPhotodisintegration_VarIRB ) {
    _fPhotodisintegration_VarIRBFlag      = 1 ;

    //get position in photon field (for table creation)
    TiXmlElement* lPD_VarIRB_X_Pos = lpXmlCheckPhotodisintegration_VarIRB->FirstChildElement("Pos_X_Mpc");
    if(lPD_VarIRB_X_Pos){
      Positon.setX(atof(GetAttribute(lPD_VarIRB_X_Pos,"value").c_str())*Mpc);
    }else{ Positon.setX(0.*Mpc); }

    TiXmlElement* lPD_VarIRB_Y_Pos = lpXmlCheckPhotodisintegration_VarIRB->FirstChildElement("Pos_Y_Mpc");
    if(lPD_VarIRB_Y_Pos){
      Positon.setY(atof(GetAttribute(lPD_VarIRB_Y_Pos,"value").c_str())*Mpc);
    }else{ Positon.setY(0.*Mpc); }

    TiXmlElement* lPD_VarIRB_Z_Pos = lpXmlCheckPhotodisintegration_VarIRB->FirstChildElement("Pos_Z_Mpc");
    if(lPD_VarIRB_Z_Pos){
      Positon.setZ(atof(GetAttribute(lPD_VarIRB_Z_Pos,"value").c_str())*Mpc);
    }else{ Positon.setZ(0.*Mpc); }
    
    std::cout<<"Position in Mpc is (x,y,z)=("
	     <<Positon.getX()/Mpc<<","
	     <<Positon.getY()/Mpc<<","
	     <<Positon.getZ()/Mpc<<")"<<std::endl;
   
  } else _fPhotodisintegration_VarIRBFlag = 0 ;

  TiXmlElement* lpXmlCheckDecayFlag = lpXmlInter->FirstChildElement("NoDecay") ;
  if ( !lpXmlCheckDecayFlag ) {
    _fDecayFlag = 1 ;
  } else _fDecayFlag = 0 ;
  TiXmlElement* lpXmlCheckCutcascade = lpXmlInter->FirstChildElement("CutCascade_MagField") ;
  if ( !lpXmlCheckCutcascade ) {
    _fCutcascadeFlag = 0 ;
  } else {
    const char* lDum = lpXmlCheckCutcascade->Attribute("value",&_fCutcascadeFlag);
    if ( lDum == NULL ) _fCutcascadeFlag = 1 ;
  }
  TiXmlElement* lpXmlCheckPP = lpXmlInter->FirstChildElement("PPProd") ;   // LM
  if ( lpXmlCheckPP ) {
    _fPProdFlag = 1 ;
  } else _fPProdFlag = 0 ;  

  

  //Check for create table flag
  string lGetCreateTableType;
  TiXmlElement* lpXmlGetCreateTable = lpXmlInter->FirstChildElement("CreatePDTable") ;
  if(lpXmlGetCreateTable){
    lGetCreateTableType = GetAttribute(lpXmlGetCreateTable,"type") ;
    if(lGetCreateTableType=="MFPTable" || lGetCreateTableType=="AveragedCrossSectionTable"){
      _fCreatePDTableFlag=1;
    }else throw TXmlErr("You've selected CreatePDTable: Define the type of table: MFPTable or AveragedCrossSectionTable?") ;
  }else _fCreatePDTableFlag=0;

  //Check for reduce table flag
  string lGetReduceTableType;
  TiXmlElement* lpXmlGetReduceTable = lpXmlInter->FirstChildElement("ReducePDTable") ;
  if(lpXmlGetReduceTable){
    lGetReduceTableType = GetAttribute(lpXmlGetReduceTable,"type") ;
    if(lGetReduceTableType=="MFPTable" || lGetReduceTableType=="AveragedCrossSectionTable"){
      _fReducePDTableFlag=1;
    }else throw TXmlErr("You've selected ReducePDTable: Define the type of table: MFPTable or AveragedCrossSectionTable?") ;
  }else _fReducePDTableFlag=0;

  int NIntegrationCalls;
  TiXmlElement* lpXmlNIntegrationCalls = lpXmlInter->FirstChildElement("NIntegrationCalls") ;
  if (lpXmlNIntegrationCalls){ 
    const char* lNIntegrationCalls = lpXmlNIntegrationCalls->Attribute("value",&NIntegrationCalls) ;
  }else{
    NIntegrationCalls=30;
  }

  //photodisintegration settings.
  //default path
  std::ostringstream TALYSDirBuffer;
  TALYSDirBuffer<<DEFAULT_TALYSTabMeanFreePath_DIR;
   TalysDirectory =  TALYSDirBuffer.str();
  
  //here is the place to define the photonfields (additonal one to CMB/IRO) that should be taken into account for PD.
  bool DefaultPDTablesSet=false;
  TiXmlElement* lpXmlCheckPDisTable1 = lpXmlInter->FirstChildElement("PhotoDisMFPTabulated");
  
  if(_fPhotodisintegration_VarIRBFlag==1 && lpXmlCheckPDisTable1){
    throw TXmlErr("You can not choose the variable IRB in combination with the MFP tables.") ;
  }
  
  //case: PD-MFP tables
  if(lpXmlCheckPDisTable1){
    
    std::cout<<"Mean free path tables for photodisintegration are choosen."<<std::endl;
    std::cout<<"You chose to use the tabulated mean free path tables for photo-disintegration."
	     <<" Note, that this assumes constant photon fields"
	     <<std::endl;
    if(_fRedshiftFlag == 1){
      std::cout<<"Note: redshift is switched on, but you chose to use the mean free path tables."
	       <<" Therefore, the MFP for PD in IR background will be scaled as a CMB-like photon field."
	       <<std::endl;
    }

    if(_fCMBPhotodisintegrationFlag  == 1 && _fIRPhotodisintegrationFlag == 0){
      //CMB for Photodisintegration
      std::cout<<"\t-> CMB mean free path tables seletced."<<std::endl;
      std::ostringstream buffer;
      buffer<<DEFAULT_TALYSTabMeanFreePath_DIR;
      TalysDirectory = buffer.str();
      DefaultPDTablesSet=true;
      _fTabulatedTALYSY_Vec.push_back(TabulatedTALYSMeanFreePath_CMB::GetInstance());
    }else if(_fIRPhotodisintegrationFlag == 1 && _fCMBPhotodisintegrationFlag  == 0){  
      //IRO for Photodisintegration
      std::cout<<"\t-> IRB mean free path tables seletced."<<std::endl;
      std::ostringstream buffer;
      buffer<<DEFAULT_TALYSTabMeanFreePath_DIR;
      TalysDirectory = buffer.str();
      DefaultPDTablesSet=true;
      _fTabulatedTALYSY_Vec.push_back(TabulatedTALYSMeanFreePath_IRB::GetInstance());
    }else if(_fIRPhotodisintegrationFlag == 1 && _fCMBPhotodisintegrationFlag  == 1){   
      //CMB+IRO for Photodisintegration
      std::cout<<"\t-> CMB and IRB mean free path tables seletced"<<std::endl;
      std::ostringstream buffer;      
      buffer<<DEFAULT_TALYSTabMeanFreePath_DIR;
      TalysDirectory = buffer.str();
      DefaultPDTablesSet=true;
      std::cout<<"\t1. TalysDirectory="<<TalysDirectory<<std::endl;
      _fTabulatedTALYSY_Vec.push_back(TabulatedTALYSMeanFreePath_CMB::GetInstance());
      buffer.str("");
      TalysDirectory.clear(); 
      buffer<<DEFAULT_TALYSTabMeanFreePath_DIR;
      TalysDirectory = buffer.str();
      std::cout<<"\t2. TalysDirectory="<<TalysDirectory<<std::endl;
      _fTabulatedTALYSY_Vec.push_back(TabulatedTALYSMeanFreePath_IRB::GetInstance());
    }else throw TXmlErr("No mean free path table directory specified in sophiainteractions.cc.") ;

    if(_fPhotonBackgrounds_Vec.size()>0){
      std::cout<<"You want to use a photonfield and the mean free path table?"
	       <<" Note, that CMB/IRO were already taken into account while creating "
	       <<"the standard mean free path tables."
	       <<std::endl;
    }
    _fTALYSMeanFreePathCalc = new TALYSMeanFreePathTabulated();
    _fTALYSMeanFreePathCalc->SetNIntegratorCalls(NIntegrationCalls);
  }
  
  //case: PD-XS tables
  double epsilon0 = 0., epsilon_min = 0., epsilon_max = 0.;
  lpXmlCheckPDisTable1 = lpXmlInter->FirstChildElement("PhotoDisAvrgdXS");
  if(lpXmlCheckPDisTable1){
    DefaultPDTablesSet=true;
    std::cout<<"Avrg xs tables choosen."<<std::endl;

    TiXmlElement* lPD_VarIRB_EpsMin = lpXmlCheckPDisTable1->FirstChildElement("epsilon_min");
    if(lPD_VarIRB_EpsMin){
     epsilon_min = atof(GetAttribute(lPD_VarIRB_EpsMin,"value").c_str());
    }else{ epsilon_min = 4.e-19; }
    
    TiXmlElement* lPD_VarIRB_EpsMax = lpXmlCheckPDisTable1->FirstChildElement("epsilon_max");
    if(lPD_VarIRB_EpsMax){
     epsilon_max = atof(GetAttribute(lPD_VarIRB_EpsMax,"value").c_str());
    }else{ epsilon_max = 12.4e-9; }

    std::cout<<"Adjustment for the numerical calculation of PD mean free path:\n"
	     <<"epsilon_min="<<epsilon_min
	     <<"epsilon_max="<<epsilon_max
	     <<std::endl;    

    if(_fRedshiftFlag == 0){
      std::cout<<"Warning: redshift is switched off, but you chose to use the averaged cross sections tables."
	       <<" Have you considered using the mean free path tables for faster simulations?"
	       <<std::endl;
    }
    
    if(_fCMBPhotodisintegrationFlag  == 0 && _fIRPhotodisintegrationFlag == 0){
      std::cout<<"Warning: You choose <PhotoDisMFPTabulated/> for Photodisintegration but switched of CMB/IR. "<<std::endl;
      std::cout<<"Did you choose a different, additional photonfield?"<<std::endl;
    }
    //TODO NN: To be changed: the photonbackgrounds should be take from the Universe* instead!
    //CMB for Photodisintegration
    if(_fCMBPhotodisintegrationFlag  == 1){
      std::cout<<"\t->Add CMB to list of activated photonfields."<<std::endl;
      _fPhotonBackgrounds_Vec.push_back(new CMB());
    }
    //IRO for Photodisintegration
    if(_fIRPhotodisintegrationFlag == 1){  
      std::cout<<"\t->Add IR to list of activated photonfields."<<std::endl;
      _fPhotonBackgrounds_Vec.push_back(new IR());
    }   
    //variable IRO for Photodisintegration
    if(_fPhotodisintegration_VarIRBFlag == 1){  
      std::cout<<"\t->Add variable IR to list of activated photonfields.aFileName:"<<aFileName<<std::endl;

      TiXmlElement* lpXmlUniverse = XmlExtract().GetElement("Environment") ;
      string lUnivType = lpXmlUniverse->Attribute("type") ;
      if (lUnivType == "One Dimension") {
	_fPhotonBackgrounds_Vec.push_back(new TVariableInfrared(aFileName, true));
      }  else if (lUnivType == "LSS")
	_fPhotonBackgrounds_Vec.push_back(new TVariableInfrared(aFileName, false));
      else
	TXmlErr("Wrong Environment type during IRB initialization for photo disintegration in sophiacinteractions.cc") ;
      
    }
    
      _fTabulatedTALYSY_Vec.push_back(TabulatedTALYSAveragedCrossSection::GetInstance());
      if(_fCreatePDTableFlag==0){
	_fTALYSMeanFreePathCalc = new TALYSMeanFreePathAvrgd(epsilon0, epsilon_min, epsilon_max);
	_fTALYSMeanFreePathCalc->SetNIntegratorCalls(NIntegrationCalls);
      }
  }
  
  if(DefaultPDTablesSet==false && _fCreatePDTableFlag==0){
    std::cout<<"Set default table for photodisintegration (MFP-tables in CMB and IRB)."<<std::endl;
    std::ostringstream buffer;      
    buffer<<DEFAULT_TALYSTabMeanFreePath_DIR;
    TalysDirectory = buffer.str();
    _fTabulatedTALYSY_Vec.push_back(TabulatedTALYSMeanFreePath_CMB::GetInstance());
    buffer.str("");
    TalysDirectory.clear(); 
    buffer<<DEFAULT_TALYSTabMeanFreePath_DIR;
    TalysDirectory = buffer.str();
    _fTabulatedTALYSY_Vec.push_back(TabulatedTALYSMeanFreePath_IRB::GetInstance());
    _fTALYSMeanFreePathCalc = new TALYSMeanFreePathTabulated();
    _fTALYSMeanFreePathCalc->SetNIntegratorCalls(NIntegrationCalls);
  }

  ////////////////////////////////////////////////////////////////
  //Creation of the PD tables 
  //1. create MFP Table
  if (_fCreatePDTableFlag==1 && lGetCreateTableType=="MFPTable" ) { 
      
      TiXmlElement* folder = lpXmlInter->FirstChildElement("CreationDirectory") ;
      TiXmlNode* folderName = folder->FirstChild() ;
      std::cout<<"Try to create the TALYS mean free path tables in directory: "
  	       <<folderName->Value()<<std::endl;
      
      //NIntegratorCalls
      TiXmlElement* NIntegratorCalls = lpXmlInter->FirstChildElement("NIntegratorCalls") ;
      int NIntegratorCallsVal;
      if(NIntegratorCalls){
  	NIntegratorCalls->Attribute("value", &NIntegratorCallsVal) ;
      }else{
  	NIntegratorCallsVal=30;
      }
      std::cout<<"NIntegratorCallsVal="<<NIntegratorCallsVal<<std::endl;
	    
      //EMinMFPTable
      TiXmlElement* EMinMFPTable = lpXmlInter->FirstChildElement("EMinMFPTable") ;
      double EMinMFPTableVal;
      if(EMinMFPTable){
  	EMinMFPTable->Attribute("value", &EMinMFPTableVal) ;
      }else{
  	EMinMFPTableVal=6;
      }
      std::cout<<"EMinMFPTableVal="<<EMinMFPTableVal<<std::endl;

      //EMaxMFPTable
      TiXmlElement* EMaxMFPTable = lpXmlInter->FirstChildElement("EMaxMFPTable") ;
      double EMaxMFPTableVal;
      if(EMaxMFPTable){
  	EMaxMFPTable->Attribute("value", &EMaxMFPTableVal) ;
      }else{
  	EMaxMFPTableVal=14;
      }
      std::cout<<"EMaxMFPTableVal="<<EMaxMFPTableVal<<std::endl;

      //NBinsMFPTable
      TiXmlElement* NBinsMFPTable = lpXmlInter->FirstChildElement("NBinsMFPTable") ;
      int NBinsMFPTableVal;
      if(NBinsMFPTable){
  	NBinsMFPTable->Attribute("value", &NBinsMFPTableVal) ;
      }else{
  	NBinsMFPTableVal=200;
      }
      std::cout<<"NBinsMFPTableVal="<<NBinsMFPTableVal<<std::endl;
      
                  
      if(_fCMBPhotodisintegrationFlag  == 1 && _fIRPhotodisintegrationFlag == 0){
      TabulatedTALYSMeanFreePath_CMB::GetInstance()->CreateTables(_fPhotonBackgrounds_Vec,
								  _fTabulatedTALYSY_Vec,
								  _fTabulatedTALYSY_Vec[0]->GetPD_TabDataType(),
								  Positon,
								  epsilon0, 
								  epsilon_min, 
								  epsilon_max,
								  folderName->Value(),
								  NIntegratorCallsVal,
								  EMinMFPTableVal,
								  EMaxMFPTableVal,
								  NBinsMFPTableVal);
      }
      if(_fCMBPhotodisintegrationFlag  == 0 && _fIRPhotodisintegrationFlag == 0 && _fPhotodisintegration_VarIRBFlag == 1){
      TabulatedTALYSMeanFreePath_IRB::GetInstance()->CreateTables(_fPhotonBackgrounds_Vec,
								  _fTabulatedTALYSY_Vec,
								  _fTabulatedTALYSY_Vec[0]->GetPD_TabDataType(),
								  Positon,
								  epsilon0, 
								  epsilon_min, 
								  epsilon_max,
								  folderName->Value(),
								  NIntegratorCallsVal,
								  EMinMFPTableVal,
								  EMaxMFPTableVal,
								  NBinsMFPTableVal);
      }
      exit(0);
  }

  //2. create averaged cross section table
  if (_fCreatePDTableFlag==1 && lGetCreateTableType=="AveragedCrossSectionTable" ) { 
      TiXmlElement* folder = lpXmlInter->FirstChildElement("CreationDirectory") ;
      TiXmlNode* folderName = folder->FirstChild() ;
      std::cout<<"Try to create the tabulated TALYS averaged cross section tables in directory: "
	       <<folderName->Value()<<std::endl;
      TabulatedTALYSAveragedCrossSection::GetInstance()->CreateTables(folderName->Value());
      exit(0);
  }

  ////////////////////////////////////////////////////////////////
  //Reduce the PD tables 
  //1. reduce MFP table
  if (_fReducePDTableFlag==1 && lGetReduceTableType=="MFPTable" ) {  

    TiXmlElement* folder = lpXmlInter->FirstChildElement("CreationDirectory") ;
    TiXmlNode* folderName = folder->FirstChild() ;
    
    TiXmlElement* AccuracyTag = 
      lpXmlInter->FirstChildElement("Accuracy_percent") ;
    double Accuracy;
    AccuracyTag->Attribute("value",&Accuracy) ;
    
    std::cout<<"Try to create the reduced tabulated TALYS MFP in directory: "
	     <<folderName->Value()<<std::endl;
    std::cout<<"Accuracy="<<Accuracy<<std::endl;
    if(_fCMBPhotodisintegrationFlag  == 1 && _fIRPhotodisintegrationFlag == 0){
    TabulatedTALYSMeanFreePath_CMB::GetInstance()->ReduceTables(folderName->Value(), 
								Accuracy);
    }
    if(_fCMBPhotodisintegrationFlag  == 0 && _fIRPhotodisintegrationFlag == 1){
      TabulatedTALYSMeanFreePath_IRB::GetInstance()->ReduceTables(folderName->Value(), 
								  Accuracy);
    }
    exit(0);
  }

  //2. reduce cross section table
  if (_fReducePDTableFlag==1 && lGetReduceTableType=="AveragedCrossSectionTable" ) { 
    TiXmlElement* folder = lpXmlInter->FirstChildElement("CreationDirectory") ;
    TiXmlNode* folderName = folder->FirstChild() ;
    
    TiXmlElement* AccuracyTag =  lpXmlInter->FirstChildElement("Accuracy_percent") ;
    double Accuracy;
    AccuracyTag->Attribute("value",&Accuracy) ;
    std::cout<<"Try to create the reduced tabulated TALYS averaged cross section tables in directory: "
	     <<folderName->Value()<<std::endl;
    
    
    std::cout<<"Accuracy="<<Accuracy<<std::endl;
    TabulatedTALYSAveragedCrossSection::GetInstance()->ReduceTables(folderName->Value(), 
								    Accuracy);
    exit(0);
  }

  
  if (!_fPionProdFlag && !_fPairProdFlag && !_fIRPionProdFlag && !_fPProdFlag && !_fPhotodisintegrationFlag) 
    std::cout<<"Warning:Sophia interactions but all interactions switched off."<<std::endl; 

    //throw TXmlErr("Sophia interactions but all interactions switched off.") ;
  //if (!_fPionProdFlag && !_fPProdFlag && (_fSecPhotonFlag || _fSecNuFlag))
  //throw TXmlErr("Secondary photons/nu but no pionproduction.") ;
  //if (!_fPairProdFlag && _fSecPairProdPhotonFlag)
  //  throw TXmlErr("Secondary photons but no pair production.") ;

  // Directory for shower tables
  if (_fSecPhotonFlag || _fSecPairProdPhotonFlag) {
    TiXmlElement* lpXmlShowerDir = lpXmlInter->FirstChildElement("ShowerTableDirectory") ;
    //std::cout<< "lpXmlShowerDir " << lpXmlShowerDir << std::endl;
    if ( ! lpXmlShowerDir ) {
//      throw TXmlErr("Incorrect photon shower table directory");
      _fShowerTableDir = DEFAULT_SHOWER_DIR ;
    }
    else {
      TiXmlNode* lpXmlShowerDirName = lpXmlShowerDir->FirstChild() ;
      if ( (! lpXmlShowerDirName) || (! lpXmlShowerDirName->ToText()) )
        throw TXmlErr("Incorrect photon shower table directory name");
      _fShowerTableDir = lpXmlShowerDirName->Value() ;
    }
  }
 
  

  // Reading interaction tables
  char cha[256] ;
  double a, b;

  //string lIntFile = lIntDir+"pairprodrate_extended.txt" ; 
  string lIntFile = lIntDir+"pair_rate_cmbir.table";

  // Added Aug 07: if the user asks no pion prod on IR we assume the user doesn't
  // want pair prod on IR either and use a dedicated pair prod table.
  if ( !_fIRPionProdFlag) lIntFile = lIntDir+"pair_rate_cmb.table";
  //  if ( !_fIRPionProdFlag) lIntFile = lIntDir+"pair_rate.table";pair_rate_cmbir.table
  fstream totppp(lIntFile.c_str(), ios::in) ;
  if (! totppp.is_open()) throw TCrpErr("Error opening Sophia interaction file : " + lIntFile);
  totppp.getline(cha,255) ;
  totppp.getline(cha,255) ;
  totppp.getline(cha,255) ;
  for (int i=0; i<70; i++) {
    totppp >> a >> b ;
    _fE_part.push_back(a / 1.e18) ; // f77 units
    _fL_pair.push_back(b/a) ; // f77 units
  }
  if (!totppp.good()) throw TCrpErr("Error reading pairprodrate_extended[_cmb].txt") ;
  totppp.close() ;
  _fdEtabbin = log10( _fE_part.at(69)/_fE_part.at(0) )/69 ;

  lIntFile = lIntDir+"pionprodrate_p.txt" ;
  fstream pionproton(lIntFile.c_str(), ios::in) ;
  if (! pionproton.is_open()) {
    throw TCrpErr("Error opening Sophia interaction file : " + lIntFile);
  }
  pionproton.getline(cha,100) ;
  pionproton.getline(cha,100) ;
  pionproton.getline(cha,100) ;
  for (int i=0; i<100; i++) {
    pionproton >> a >> b ;
    _fE_pionprod.push_back(a * EeV) ;
    _fLossRateProton.push_back(b / Mpc) ; // table is in Mpc^(-1)
  }
  if (!pionproton.good()) throw TCrpErr("Error reading pionprodrate_p.txt") ;
  pionproton.close() ;
  _fdEtabPion = log10( _fE_pionprod.at(99)/_fE_pionprod.at(0) )/99 ;

  lIntFile = lIntDir+"pionprodrate_n.txt" ;
  fstream pionneutron(lIntFile.c_str(), ios::in) ;
  if (! pionneutron.is_open()) {
    throw TCrpErr("Error opening Sophia interaction file : " + lIntFile);
  }
  pionneutron.getline(cha,100) ;
  pionneutron.getline(cha,100) ;
  pionneutron.getline(cha,100) ;
  for (int i=0; i<100; i++) {
    pionneutron >> a >> b ;
    _fLossRateNeutron.push_back(b / Mpc) ; // table is in Mpc^(-1)
  }
  if (!pionneutron.good()) throw TCrpErr("Error reading pionprodrate_n.txt") ;
  pionneutron.close() ;

  lIntFile = lIntDir+"pionprodrate_p_ir_z0.txt" ;
  fstream irpionproton(lIntFile.c_str(), ios::in) ;
  if (! irpionproton.is_open()) {
    throw TCrpErr("Error opening Sophia interaction file : " + lIntFile);
  }
  irpionproton.getline(cha,100) ;
  irpionproton.getline(cha,100) ;
  irpionproton.getline(cha,100) ;
  for (int i=0; i<150; i++) {
    irpionproton >> a >> b ;    
    _fE_IRpionprod.push_back(a * EeV);
    _fIRLossRateProton.push_back(b / Mpc) ; // table is in Mpc^(-1)
    
    //    std::cout << a << " failbit " << irpionproton.fail() << " badbit " << irpionproton.bad() << " eof " << irpionproton.eof() << std::endl;
  }
  //_fdEtabIRPion  = log10( _fE_IRpionprod.at(99)/_fE_IRpionprod.at(0) )/99. ;
  //std::cout<< _fdEtabIRPion << std::endl;
  _fdEtabIRPion  = log10( _fE_IRpionprod.back()/_fE_IRpionprod.front() )/double(_fE_IRpionprod.size()-1);
  //std::cout<< _fdEtabIRPion << std::endl;
  if (irpionproton.fail()){
    std::cout << "failbit " << irpionproton.fail() << " badbit " << irpionproton.bad() << std::endl;
 throw TCrpErr("Error reading pionprodrate_p_ir_z0.txt") ;
  }
  irpionproton.close() ;
  

  lIntFile = lIntDir+"pionprodrate_n_ir_z0.txt" ;
  fstream irpionneutron(lIntFile.c_str(), ios::in) ;
  if (! irpionneutron.is_open()) {
    throw TCrpErr("Error opening Sophia interaction file : " + lIntFile);
  }
  irpionneutron.getline(cha,100) ;
  irpionneutron.getline(cha,100) ;
  irpionneutron.getline(cha,100) ;
  for (int i=0; i<150; i++) {
    //irpionneutron >> a >> b ;
    //while(!irpionneutron.eof()){	
    _fIRLossRateNeutron.push_back(b / Mpc ) ; // table is in Mpc^(-1)
    irpionneutron >> a >> b ;
  }
  if (irpionneutron.fail()) throw TCrpErr("Error reading pionprodrate_n_ir_z0.txt") ;
  irpionneutron.close() ;

  // Neutron decay table : energy distribution of outcoming electron
  double lQ = neutron_mass_c2 - proton_mass_c2 ;
  int lSize = 100 ;
  double* lEnergies = new double[lSize] ;
  double* lDensity = new double[lSize] ;
  for (int i=0; i<lSize; i++) {
    lEnergies[i] = electron_mass_c2 + i*(lQ-electron_mass_c2)/(lSize-1) ;
    lDensity[i] = lEnergies[i]*
      sqrt(max(lEnergies[i]*lEnergies[i]-electron_mass_c2*electron_mass_c2,0.))*
      (lQ-lEnergies[i])*(lQ-lEnergies[i]) ;
  }
  _fpRandDistriNeutronDecay = new RandGeneral(lDensity,lSize,0) ;
  delete[] lDensity ;
  delete[] lEnergies ;

//   //  Pair production table : in the same directory as fShowerTableDir.
//   string lPairProdTable1 = _fShowerTableDir+"/pairprod_table.txt";
//   fstream pptable1(lPairProdTable1.c_str(),ios::in);
//   double _lPairProdTable1[60][170];
//   if (! pptable1.is_open()) throw TCrpErr("Error opening pair production table : "+lPairProdTable1);
//   //char cha[256]; double a;
//   pptable1.getline(cha,200);
//   pptable1.getline(cha,200);
//   pptable1.getline(cha,200);
//   for (int i=0; i<60; i++) {
//     for (int j=0; j<170; j++) {
//       pptable1 >> a;
//       _lPairProdTable1[i][j]=a;
//     }
//   }
//   if (!pptable1.good()) throw TCrpErr("Error reading pairprod_table.txt") ;
//   pptable1.close() ;
//   double lEpmin1=pow(10.,17);
//   double lEpmax1=pow(10.,23);
//   double _lNucleonEnergy_PP1[60];
//   for (int i=0; i<60; i++) _lNucleonEnergy_PP1[i]=lEpmin1*pow(lEpmax1/lEpmin1,double(i)/double(60));

//   {
//     fstream outPPT("PPTold.txt", ios::out);
//     for (int i=0; i<60; i++) {
//       for (int j=0; j<170; j++) {
// 	outPPT << _lNucleonEnergy_PP1[i] << '\t' << _lPairProdTable1[i][j] << std::endl;
//       }
//     }
//     outPPT.close();
//   }
  // Pair production table from Ricard
  // It is in the shower dir
  
  if (_fSecPhotonFlag || _fSecPairProdPhotonFlag) {
    string lPairProdTable = _fShowerTableDir+"pair_spectrum_cmb.table";
    if( _fIRPionProdFlag) lPairProdTable = _fShowerTableDir+"pair_spectrum_cmbir.table";
    fstream pptable(lPairProdTable.c_str(),ios::in);
    if (! pptable.is_open()) throw TCrpErr("Error opening pair production table [Ricard]: "+lPairProdTable);

    //  char cha[256]; double a;
  //pptable.getline(cha,200);
  //pptable.getline(cha,200);
  //pptable.getline(cha,200);

    double lEp, lEe;
    //std::cout << "conversion constant : " << c_light /centimeter * second << std::endl;
    for (int i=0; i<NBINS_PAIRPROD; i++) {
      for (int j=0; j<170; j++) {
	pptable >> lEp >> lEe >> a;
	_fPairProdTable[i][j]= lEe * lEe * a * c_light /centimeter * second;
      }
      _fNucleonEnergy_PP[i]=lEp;
    }
    if (!pptable.good()) throw TCrpErr("Error reading pairprod_table.txt") ;
    pptable.close() ;
  }
  //  double lEpmin=_fNucleonEnergy_PP[0];
  //  double lEpmax=_fNucleonEnergy_PP[69];
  
//    {
//      fstream outPPT("PPTnew.txt", ios::out);
//      for (int i=0; i<70; i++) {
//        for (int j=0; j<170; j++) {
//  	outPPT << _fNucleonEnergy_PP[i] << '\t' <<_fPairProdTable[i][j] << std::endl;
//        }
//      }
//      outPPT.close();
//    }
  // LM : Starts proton proton
  if (_fPProdFlag) {
    lIntFile = lIntDir+"cross_pp.dat" ;
    fstream ppproton(lIntFile.c_str(), ios::in) ;
    if (!ppproton.is_open()) {
      throw TCrpErr("Error opening Sibyll proton interaction file : " + lIntFile);
    }
    ppproton.getline(cha,100) ;
    ppproton.getline(cha,100) ;
    ppproton.getline(cha,100) ;
 
    if (!pionproton.good()) throw TCrpErr("Error reading cross_pp.dat") ;
    while (ppproton.good()) {
    //   for (int i=0; i<100; i++) {
      ppproton >> a >> b ;
      _fE_PP.push_back(a * GeV) ;
      _fXSecPPProton.push_back(b * millibarn ) ; // table is in mb
    }
    _fE_PP.pop_back();
    _fXSecPPProton.pop_back();
    //    for (int i = 0; i < _fE_PP.size(); i++) cerr << _fE_PP[i] << " " << _fXSecPPProton[i] << endl;
    ppproton.close() ;
     
  }
  // LM : Ends proton proton
  double c;
  lIntFile = lIntDir+"sigma_pg_p.txt" ;
  fstream irpptable(lIntFile.c_str(), ios::in);
  if (! irpptable.is_open()) {
    throw TCrpErr("Error opening IR integral proton file : " + lIntFile);
  }
  for (int i=0; i<501; i++) {
    for (int j=0; j<1000; j++) {
      irpptable >> a >> b >> c;
      if (j==0) _fE_IR.push_back(a);
      if (i==0) _fS_IR.push_back(b);
      _fIRXsecInt_Proton.push_back(c) ;
    }
  }
  if (!irpptable.good()) throw TCrpErr("Error reading " + lIntFile) ;
  irpptable.close() ;
  _fdEtabXsecIR = log10(_fE_IR[1]/_fE_IR[0]);
  _fdStabXsecIR = _fS_IR[1]-_fS_IR[0];

  lIntFile = lIntDir+"sigma_pg_n.txt" ;
  fstream irpntable(lIntFile.c_str(), ios::in);
  if (! irpntable.is_open()) {
    throw TCrpErr("Error opening IR integral proton file : " + lIntFile);
  }
  for (int i=0; i<251; i++) {
    for (int j=0; j<1000; j++) {
      irpntable >> a >> b >> c;
      if (j==0) _fE_IR.push_back(a);
      if (i==0) _fS_IR.push_back(b);
      _fIRXsecInt_Neutron.push_back(c) ;
    }
  }
  if (!irpntable.good()) throw TCrpErr("Error reading " + lIntFile) ;
  irpntable.close() ;

}

double TSophiaInteractions::XSecPPProton(double en) {
  if (en > _fE_PP.back()) return _fXSecPPProton.back();

  int j = floor(log10(en/_fE_PP[0])/log10(_fE_PP[1]/_fE_PP[0]));
  int i = j+1;

  //  cerr << en << " " << _fE_PP[j] << " " << _fE_PP[i] << endl;

  //  while (_fE_PP[i] < en) i++;
  //  int j = i-1;

  double r1 = (en-_fE_PP[j])/(_fE_PP[i]-_fE_PP[j]);
  double result = _fXSecPPProton[j] * (1.0-r1) + r1 * _fXSecPPProton[i]; 
  return result;
}

TSophiaInteractions::~TSophiaInteractions() {

  for(int i=0; i<_fPhotonBackgrounds_Vec.size(); i++) delete _fPhotonBackgrounds_Vec[i];
  delete  _fTALYSMeanFreePathCalc;
  
  delete _fpRandDistriNeutronDecay ;
  //delete _fpNucleusDB;
// _fpEvolModel not deleted because it is managed by the Universe class.
}

double TSophiaInteractions::PairProdSpec(double aEp_eV, int aBin_Ee) const {
  double lInterpolValue=0;
  if (aEp_eV < _fNucleonEnergy_PP[0] || aEp_eV > _fNucleonEnergy_PP[NBINS_PAIRPROD-1]) lInterpolValue=0;
  //else if (aEp_eV > _fNucleonEnergy_PP[NBINS_PAIRPROD-1]) lInterpolValue=0;
  else {
//    int lBin_Ep=0;
      int lBin_Ep = upper_bound(_fNucleonEnergy_PP,_fNucleonEnergy_PP+NBINS_PAIRPROD,aEp_eV)-_fNucleonEnergy_PP;
      /*
      cout << "Begin" << endl;
      cout << "lBin_Ep new " << lBin_Ep << endl;
      lBin_Ep = 0;
    while (aEp_eV > _fNucleonEnergy_PP[lBin_Ep]) {
      lBin_Ep++;
      if (lBin_Ep == NBINS_PAIRPROD-2) break;
    }
      cout << "lBin_Ep old " << lBin_Ep << endl;
       */
    double alpha=(aEp_eV-_fNucleonEnergy_PP[lBin_Ep])/(_fNucleonEnergy_PP[lBin_Ep+1]-_fNucleonEnergy_PP[lBin_Ep]);
    lInterpolValue = (1-alpha)*_fPairProdTable[lBin_Ep][aBin_Ee] + alpha*_fPairProdTable[lBin_Ep+1][aBin_Ee];
  }
  
  return lInterpolValue;
}

//This function returns the Id of an exclusive channel - keep in mind that we already know that there was a reaction. 
unsigned long int TSophiaInteractions::GetExclusivePDChannelDirectly(unsigned long int InitNucleusId,
								     TVector3D Position,
								     double Energy,
								     double NucleusMass,
								     double z,
								     double TimeStep){
  while(1){// Loop until a channel is found
  //Now the selection of the exclusive Channel
  double RandomNumber = RandFlat::shoot();

  double PDAllChannelsMeanFreePath=1./this->GetPDTimeStep(InitNucleusId,
						       Position,
						       Energy,  //already redshifted*/1000
						       NucleusMass, //already *c*c/1000 
						       z);
  PDAllChannelsMeanFreePath/=Mpc * inv_c_light;
  double ActualProbability=0.;
  for(int j=0; j<_fTabulatedTALYSY_Vec.size(); j++){
    
    std::vector<unsigned long int> ListOfExclChannels = _fTabulatedTALYSY_Vec[j]->GetListOfAllExclusiveChannels(InitNucleusId);
    //ListOfExclChannels = _fTabulatedTALYSY_Vec[j]->GetListOfAllExclusiveChannels(InitNucleusId);
    
    //for(int i=0; i<ListOfExclChannels.size(); i++){
      for (std::vector<unsigned long int>::iterator it = ListOfExclChannels.begin(); it != ListOfExclChannels.end(); ++it) {
      //ignore nuclei excitation
      //if(ListOfExclChannels[i]==0) continue;
          if (!(*it)) continue;
      double FreeMeanPath = _fTALYSMeanFreePathCalc->Total_rate_log(_fPhotonBackgrounds_Vec,
								    _fTabulatedTALYSY_Vec,
								    _fTabulatedTALYSY_Vec[j]->GetPD_TabDataType(),
								    Position,
								    InitNucleusId,
								    *it,//ListOfExclChannels[i],
								    NucleusMass,
								    Energy,  //GeV
								    0.,       //s. eq 13 from CRPropa paper
								    1); //2 means using sum over exclusive channels 
      
      
      //In case of the MFP-tables we perform a cmb-like scaling for the IRB, too. Furthermore, the IRB is switched of if z>MaxRedshift.
      double scaleFactor=pow(1+z,3);
      if(_fTabulatedTALYSY_Vec[j]->GetPD_TabDataType()==IRB_MFPTab){
	scaleFactor=_fpEvolModel->GetScalingFactor(z);
	//Note: scale factor =0 means IRB switched off <=> z>MaxRedshift.
	if(scaleFactor==0.){
	  scaleFactor=1.;
	  FreeMeanPath=std::numeric_limits<double>::min();
	}
      }
      
      FreeMeanPath=FreeMeanPath*scaleFactor;
      
      FreeMeanPath=1./FreeMeanPath;
      ActualProbability += PDAllChannelsMeanFreePath/FreeMeanPath; 
      
      if(ActualProbability > RandomNumber){
          return(*it);//(ListOfExclChannels[i]);  
      }
    }//i: all exclusive channels
  }//j: components of _fTabulatedTALYSY_Vec 
  }//Loop until a channel is found
  //  std::cout<<"RandomNumber="<<RandomNumber<<" ActualProbability="<<ActualProbability<<std::endl;
  //  std::cerr<<"Warning in in unsigned long int TSophiaInteractions::GetExclusiveChannelDirectly(...). Wasn't able to find the exclusive chanel!"<<std::endl;
  return(100000);
}

//Returns sum 1 / lambda_{i} for all exclusive channel's MFP \lambda{i} of the nucleus with Id=InitNucleusId. 
//The returned unit is time^-1.
double TSophiaInteractions::GetPDTimeStep(unsigned long int InitNucleusId,
					TVector3D Position,
					double Energy,
					double NucleusMass,
					double z){

  double TotalFreeMeanPath=0.;

  for(int j=0; j<_fTabulatedTALYSY_Vec.size(); j++){
  
    
    double FreeMeanPath = _fTALYSMeanFreePathCalc->Total_rate_log(_fPhotonBackgrounds_Vec,
								  _fTabulatedTALYSY_Vec,
								  _fTabulatedTALYSY_Vec[j]->GetPD_TabDataType(),
								  //Position_test,
								  Position,
								  InitNucleusId,
								  999, //Not needed . . . value will be ignored . . . 
								  NucleusMass,
								  Energy,  //GeV
								  0.,       //s. eq 13 from CRPropa paper
								  2); //2 means using sum over exclusive channels 
      
      
       //In case of the MFP-tables we perform a cmb-like scaling for the IRB, too. Furthermore, the IRB is switched of if z>MaxRedshift.
      double scaleFactor=pow(1+z,3);
    
      if(_fTabulatedTALYSY_Vec[j]->GetPD_TabDataType()==IRB_MFPTab){
	scaleFactor=_fpEvolModel->GetScalingFactor(z);
	//Note: scale factor =0 means IRB switched off <=> z>MaxRedshift.
	if(scaleFactor==0.){
	  scaleFactor=1.;
	  FreeMeanPath=std::numeric_limits<double>::min();
	}
      }

      FreeMeanPath=FreeMeanPath*scaleFactor;
           
      TotalFreeMeanPath += FreeMeanPath;
  }//j: components of _fTabulatedTALYSY_Vec
  
  
  double TimeStep=TotalFreeMeanPath/(Mpc * inv_c_light);
  return(TimeStep);
}

std::vector<TabulatedTALYSY*> TSophiaInteractions::GetPointerToPhotoDisTables(){
  return(_fTabulatedTALYSY_Vec);
}
