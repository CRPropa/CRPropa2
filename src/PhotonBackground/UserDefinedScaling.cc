#include <UserDefinedScaling.h>

UserDefinedScaling:: UserDefinedScaling(std::string XMLFileName): TIRBzEvolutionModel::TIRBzEvolutionModel(XMLFileName){
  //Read the MaxRedshift value from xml configuration
  cout<<"test"<<endl;
  TiXmlElement* lpXmlField = XmlExtract().GetElement("IRB_MFPScalingModel") ;
  TiXmlElement* lpXmlFile = lpXmlField->FirstChildElement("File") ;
  string lFileType = lpXmlFile->Attribute("type") ;
  cout<<"test2"<<endl;
  if (lFileType == "ASCII") {

	TiXmlNode* lpXmlFileName = lpXmlFile->FirstChild() ;
	if ( !lpXmlFileName || !(lpXmlFileName->ToText()) ) {
	  throw TXmlErr("Incorrect redshift file name");
	}
	string lFileName =lpXmlFileName->Value() ;
	fstream lScalingstream(lFileName.c_str(), ios::in) ;
	if (! lScalingstream.is_open()) {
	  throw TCrpErr("Error opening IRB scaling file : " + lFileName );
	}
	char cha[256] ;
	lScalingstream.getline(cha,100) ; // 2 lines of comment
	lScalingstream.getline(cha,100) ;
	double lDumZ, lDumScale; // convention: redshift, scaling !
	while (lScalingstream >> lDumZ >> lDumScale) {
	  cout << lDumZ << lDumScale <<endl;
	  _fRedshiftArray.push_back(lDumZ) ;
	  _fScalingArray.push_back(lDumScale) ;
	  if (!lScalingstream.good()) throw TCrpErr("Error reading IRB scaling file.") ;
	}
	lScalingstream.close() ;

  } else throw TXmlErr("Error getting file type for IRB scaling" );

}

double UserDefinedScaling::GetScalingFactor(double redshift) {

  if(redshift>_fRedshiftArray[-1]) return(0.);

  //do a linear interpolation between the specified points

  int position=0;
  while (redshift<_fRedshiftArray[position]){
	  cout<<'redshift '<< _fRedshiftArray[position]<<endl;
	  position++;
  }
  double relPositionZ=(redshift-_fRedshiftArray[position-1])/
		  (_fRedshiftArray[position]-_fRedshiftArray[position-1]);
  double relScaling= (_fScalingArray[position]-_fScalingArray[position-1])*relPositionZ;
  return _fScalingArray[position]-relScaling;

}
