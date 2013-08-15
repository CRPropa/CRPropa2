#include "TabulatedTALYSCrossSections.h"

TabulatedTALYSCrossSection* TabulatedTALYSCrossSection::finstance;

TabulatedTALYSCrossSection::TabulatedTALYSCrossSection(){ 
 fPD_TabDataType=XSTab;
 fXUnit="MeV"; 
 fXDescription="photon energy in cms";
 fYUnit="mbarn";
 fYDescription="cross section" ;
 
 std::string TalysDirectory = DEFAULT_TALYS_DIR;
 
 std::cout<<"TALYS Directory="<<TalysDirectory<<std::endl;

 //open files
 std::string IntFile = TalysDirectory+"PDTotalCross.cxs";

 IntFile = TalysDirectory+"PDExclCross.cxs";

 std::ifstream ExclusiveCrossSecData(IntFile.c_str(), std::ios::in);
 if (! ExclusiveCrossSecData.is_open()) 
    throw TCrpErr("Error in TabulatedTALYSCrossSection : Can't open : "+IntFile
		 +" Note, you have to install them by starting ./GetPDCrossSections.sh");

 IntFile = TalysDirectory+"PDExclCrossId.cxs";
 std::ifstream ExclusiveCrossSecDataId(IntFile.c_str(), std::ios::in);
 if (! ExclusiveCrossSecDataId.is_open()) 
   throw TCrpErr("Error opening TALYS Cross section data file : " + IntFile
		 +" Note, you have to install them by starting ./GetPDCrossSections.sh");
 
 IntFile = TalysDirectory+"PDInitNucleusId.cxs";
 std::ifstream TALYSInitNucleusId(IntFile.c_str(), std::ios::in);
 if (! TALYSInitNucleusId.is_open()) 
   throw TCrpErr("Error opening TALYS Cross section data file : " + IntFile
		 +" Note, you have to install them by starting ./GetPDCrossSections.sh");

 //vectors
 int Type=99;
 TALYSInitNucleusId 
   >>Type
   >>fXMin 
   >>fXMax
   >>fXNSteps;

 fXStepsize = (fXMax-fXMin) / ( (double) fXNSteps);
   
 if(Type!=0){
   throw TCrpErr("The file in "+TalysDirectory+" are no proper .cxs files!");
 }

 std::cout<<"Read in cross section data for photo-disintegration . . . this may take some time."<<std::endl;
 std::cout<<"Minimum energy [MeV]: "<<fXMin<<std::endl;
 std::cout<<"Maximum energy [MeV]: "<<fXMax<<std::endl;
 std::cout<<"Number of bins: "<<fXNSteps<<std::endl;
 
 int i=0;
 double ExclCrossSecVal;
   while (ExclusiveCrossSecData 
	  >>ExclCrossSecVal) {
     fExclusiveYData_Vec.push_back(ExclCrossSecVal) ; 
     i++;
   }
  
   std::vector<unsigned long int> ActualRow(4,0);
 
   while(ExclusiveCrossSecDataId 
	 >>ActualRow[0] 
	 >>ActualRow[1]
	 >>ActualRow[2] 
	 >>ActualRow[3]){
     fTALYSExclYId_Matrix.push_back(ActualRow);
   }
   
   std::vector<unsigned long int> ActualRow1x5(3,0);   
   while(TALYSInitNucleusId 
	 >>ActualRow1x5[0] 
	 >>ActualRow1x5[1] 
	 >>ActualRow1x5[2]
	 ){
     fTALYSInitNucleusId_Matrix.push_back(ActualRow1x5);

   }
}

TabulatedTALYSCrossSection::~TabulatedTALYSCrossSection(){
  
}

