
/**
   @file    nucleusDB.cc
   @author  Joerg Kulbartz, jkulbart@mail.desy.de
   @brief   Implementation of the TNucleusDB class. See the .h file
*/

#include "nucleusdb.h"
TNucleusDB* TNucleusDB::_instance;


TNucleusDB::TNucleusDB(void){
  /*Constructor, reads the database */ 
#ifdef DEBUG_OUTPUT
  std::cout<< "TNucleusDB constructor! "<< std::endl;
#endif

  int i, j, k;
  std::string HLTDirectory = DEFAULT_HLT_DIR;

  ifstream HLTable(( HLTDirectory +"bnl_data.dat").c_str());
  ifstream PTETable(( HLTDirectory +"PTE_bnl.dat").c_str());
  //  nucleus NDat;
  string line, another_line;
  char str[256], *decstr;
  decstr=str;
  double deltaMass, tHalf; 
  int length;
  int A, Z; 
  i=0;
  j=0;
  k=0;

  if( HLTable.fail())  throw TCrpErr("Halflife Table not found" ) ;
  if( PTETable.fail())  throw TCrpErr("PTE Table not found" ) ;
      
  HLTable.getline( str, 255);
  HLTable.getline( str, 255 );

  while(!PTETable.eof()){
    PTETable >> Z >> A >> line; 
    PTE[Z]=A; 
  }



#ifdef DEBUG_OUTPUT
  std::cout << "PTE Table" << std::endl;
  for(i=1; i<29; i++) std::cout<< i << PTE[i] << std::endl; 
#endif

  i=0;
  while(!HLTable.eof()){
    nucleus NDat;
    HLTable >> NDat.A >> NDat.Z >> NDat.N >> NDat.deltaMass >> NDat.tau >> NDat.name;
    HLTable.getline( str, 255 );//End the line read above
    HLTable.getline( str, 255 );//Reads the decay string
    
    //Parsing decay string into nucleus structure
    if(NDat.tau != 0){//Decay modes only if nuclei is unstable.
      char* decstor;
      char typeplc[256];//Allocate memory
      char* type; 
      type=typeplc;
      double br;
      decstr=(char*)str;
      decstor=decstr;
      j=0;
      while(*decstr!='\0'){
	//These are just two pointers, running along the string str, which should be formated as E:50;EP:50;
	while(*decstor != ':' && *decstor!= '\0') decstor++;// jump decstor to the ':'
	if(*decstor == '\0') throw TCrpErr("NucleusDB: Malformed decay table" );
	*decstor='\0';//As a string decstr now contains everything in front of ':'
	sscanf(decstr, "%s", type);
	decstr=++decstor;//decstr and decstor are now both directly behind the ':'
	while(*decstor != ';' && *decstor != '\0') decstor++;//jump to ';'
	if(*decstor == '\0') throw TCrpErr("NucleusDB: Malformed decay table" );
	*decstor='\0';//decstr now contains everything between ':' and ';'
	br=atof(decstr);
	NDat.decayRates.push_back(br);
	NDat.decayTypes.push_back(string(type));
	decstr=++decstor;//jump behind ';' if no other channel is present, while will abort.
      }
    }

    
    for(int j=1; j<NDat.decayRates.size(); j++) NDat.decayRates[j]+=NDat.decayRates[j-1];//Sum over the decay rates
    /*
    std::cout << NDat.A << '\t' <<  NDat.name << '\t';
    for(int j=0; j<NDat.decayRates.size(); j++) std::cout<< NDat.decayRates[j] << '\t';
    std::cout<<std::endl;
    std::cout<< "Culmulative decay rate:"<< std::endl;

    for(int j=0; j<NDat.decayRates.size(); j++) std::cout<< NDat.decayTypes[j] << '\t';
    std::cout<<std::endl;
    */
    //    if(i > 20) exit(0);
    //    i++;



      /*	
    if(NDat.decayMode != 0) {
      NDat.tau = atof(another_line.c_str());
    }else NDat.tau =0.;
      */
    NDat.tau /= log(2.);

    ret = NucleusDB.insert ( pair<int, nucleus>((NDat.A + NDat.Z * 1000), NDat ) );
    k = (NDat.A + NDat.Z * 1000);

    //std::cout<< k << " inserted" << std::endl;
    if(ret.second==false) {
      std::cout<< "Insertion failed: k=" << k << std::endl;
      throw TCrpErr("NucleusDB not able to create decay table." );
    }


    //#ifdef DEBUG_OUTPUT
    /*Complicated errorchecking lines 
    if(ret.second == false &&( NDat.decayRates[0] != NucleusDB[(NDat.A + NDat.Z * 1000)].decayRates[0] || NucleusDB[k].tau != NDat.tau || NucleusDB[k].N != NDat.N || NucleusDB[k].deltaMass != NDat.deltaMass) ){
      j++;
      cout<< " A " << NDat.A << " Z " << NDat.Z << " N " << NDat.N << " deltaMass " << NDat.deltaMass << " tHalf " << NDat.tau << " decMode Types" << NDat.decayTypes[0] << " decRates " << " DecayRates "  << NDat.decayRates[0] <<  " Name " << NDat.name << endl;
      cout<< " A " << NucleusDB[(NDat.A + NDat.Z * 1000)].A << " Z " << NucleusDB[(NDat.A + NDat.Z * 1000)].Z << " N " << NucleusDB[(NDat.A + NDat.Z * 1000)].N << " deltaMass " << NucleusDB[(NDat.A + NDat.Z * 1000)].deltaMass << " tHalf " << NucleusDB[(NDat.A + NDat.Z * 1000)].tau << " decMode " << NucleusDB[(NDat.A + NDat.Z * 1000)].decayTypes[0] << " decayRates " << NucleusDB[(NDat.A + NDat.Z * 1000)].decayRates[0] << " Name " << NucleusDB[(NDat.A + NDat.Z * 1000)].name << endl;
      } End if Errorchecking */
    //#endif
  
  }//End of while
 
}



/*Here the nucleus decay stuff starts. */ 

inline int TNucleusDB::GetNucleusKey(int A, int Z)
{
  return Z * 1000 + A;
}

inline int TNucleusDB::GetA(int aKey)
{
  it=NucleusDB.find(aKey);
  if (it != NucleusDB.end()) {
  return (*it).second.A;
  } else return aKey % 1000;
}

inline int TNucleusDB::GetZ(int aKey)
{
  it=NucleusDB.find(aKey);
  if (it != NucleusDB.end()) {
    return (*it).second.Z;
  }else return int(aKey)/int(1000);
}

inline int TNucleusDB::GetN(int aKey)
{
  it=NucleusDB.find(aKey);
  if (it != NucleusDB.end()) {
    return (*it).second.N;
  }else return (this->GetA(aKey) - int(aKey)/int(1000));
}

double TNucleusDB::GetMass(int aKey)
{
  it=NucleusDB.find(aKey);
  if (it != NucleusDB.end()) {

#ifdef DEBUG_OUTPUT
    static int flag=0;
    std::cout<< "GetMass called: aKey: "<< aKey << std::endl; 
    if (flag != aKey) {
      cout << "GetMass: aKey" << aKey << " Mass "<<( (double(aKey % 1000) * 931.5 + it->second.deltaMass - double(int(aKey / 1000))/2.) * MeV)/c_squared * c_squared << std::endl;}
    flag=aKey; 
#endif
    return  (double(aKey % 1000) * 931.5 + it->second.deltaMass - double(int(aKey / 1000)) * 511. * keV)/c_squared;}
  else{
    /*
    std::cout<< "Fallback Pass in TNucleusDB:: GetMass()"<< std::endl;
    std::cout<< "aKey :" << aKey << std::endl;
    */
    return double(aKey % 1000) * 931.5  * MeV/c_squared;
  }
}

string TNucleusDB::GetName(int aKey)
{
  it=NucleusDB.find(aKey);
  if (it != NucleusDB.end()) {
    return (*it).second.name;
  }else return "";
}

double TNucleusDB::GetTau(int aKey)
{
  it=NucleusDB.find(aKey);
  if (it != NucleusDB.end()) {
    return (*it).second.tau * second;
  }
  return 1.e-20 * second;// Nucleon dripping in this case.
  
}

string TNucleusDB::GetDecayMode(const int aKey,const double proba )
{
  /*
  static int second=0;
  static int total=0;
  static const int tstKey=16043;
*/
  it=NucleusDB.find(aKey);

  if(it!=NucleusDB.end()){
    //if(aKey==tstKey) total++;
    int k=0;
    while((it->second.decayRates[k] < proba) && (it->second.decayRates.size() < k))k++;//TODO Test this!!!
    /*
    if(k !=0 ){
      if(aKey==tstKey)second++;
      std::cout << "GetDecayMode : aKey " << aKey << std::endl;
      std::cout<< "proba " << proba << " decayRates " << it->second.decayRates[k] << std::endl;
	for(int i=0; i<=k;i++) std::cout << " decay Types " <<it->second.decayTypes[i] << " decayRates " << it->second.decayRates[i] << std::endl;
	std::cout<< "\n" << "k " << k << " second ratio : " << second/(double)total << std::endl;
    }
    */
    return it->second.decayTypes[k];//choose a decay rate
      } else {
    int Z=aKey/(int)1000;
    int A=aKey%1000;
    int N=(A-Z);
    if(A > PTE[Z]) return string("N"); //neutron dripping
    else if (A < PTE[Z]) return string("P");//proton dripping
    if(Z==0) return string("N"); //neutron dripping
    if(N==0) return string("P"); //proton dripping
    else  throw TCrpErr("Nucleus not found in GetDecayMode" );
  }
}


