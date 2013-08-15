/**
   @file    NucleusDB.h
   @author  Joerg Kulbartz, jkulbart@mail.desy.de
   @brief   Stores the properties of nuclei
*/

#ifndef _NucleusDB_H_
#define _NucleusDB_H_

#include<fstream>
#include<map>
#include<math.h>
#include"units.h"
#include"typedclass.h"
#include"crp_err.h"
#include<cstring>
#include<vector>

using namespace std;

/**
   @class TNucleusDB
   @brief Stores the properties of nuclei

*/

/* Again a singelton (compare TabulatedTalysCrossections.h)
class N
 {
 public:
    static N* instance()
    {
       if(! _instance)
       _instance = new N ();
       return _instance;
    }

    ~N();
    void xyz();

 private:
 static N* _instance;
    N();            // no instance of N can be created from the "outside" of N
                    // protected, wenn man von der Klasse noch erben möchte
    N( const N& );              //no new instance via copy constructor
    N & operator = (const N &); //no new instance via copy constructor

*/

struct nucleus{
  string name;
  int A;
  int Z; 
   int N;
  double tau;
  double deltaMass;
  //  int decayMode;
  vector<string>decayTypes;
  vector<double>decayRates;
};

class TNucleusDB : public TTypedClass {


 public:

static TNucleusDB* GetInstance()
{
  if(! _instance)  _instance = new TNucleusDB();
  return _instance;
}

  int GetA(int aKey);
  int GetZ(int aKey);
  int GetN(int aKey);
  string GetName(int aKey);
  double GetTau(int aKey);
  string GetDecayMode(int aKey, double proba);
  double GetMass(int aKey);
  int GetNucleusKey(int A, int Z);

 private: 
  /* The nuclear database */
  static TNucleusDB* _instance;
  TNucleusDB(void);
  map<int, nucleus>NucleusDB;
  map<int, nucleus>::iterator it;

  int PTE[29];
  pair<map<int, nucleus>::iterator,bool> ret;

};

#endif
