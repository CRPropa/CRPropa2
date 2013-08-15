/**
   @file   basicparam.h
   @author Tristan Beau, beau@in2p3.fr
   @brief  Class describing basic parameters for the execution of the simulation
*/

#ifndef _BASICPARAM_H_
#define _BASICPARAM_H_

#include "xmlparam.h"
#include "outputdata.h"
#include "nooutput.h"
#include "asciioutput.h"
#include "fitsoutput.h"
#include "rootoutput.h"
#include "crp_err.h"
#include "xml_err.h"

#include <stdlib.h>
#include <fstream>
#include "units.h"
#include "sysdep.h"

#include "CLHEP/Random/RandFlat.h"

using namespace std;

/**
   @class TBasicParam
   @brief Class describing basic parameters in the simulation. Those include the number of trajectories
   (primaries and secondaries), the minimum energy for nucleon propagation, the maximum time for nucleon propagation,
   the name of outputfile and record mode (full trajectories or events), a pointer to an object of the TOutputData class,
   and a random seed.
*/

class TBasicParam : public TXmlParam {
 public:
  TBasicParam(const char*) ;
  ~TBasicParam() ;

  unsigned long N() const { return _fN; } /**< Number of trajectories simulated */
  double Emin() const { return _fEmin; } /**< Minimum energy for propagation */
  double Tmax() const { return _fTmax; } /**< Maximum propagation time */

  string OutputFileName() const { return _fOutputFileName; } /**< Name of output file */
  string RecordMode() const { return _fRecordMode; } /**< Record mode (full trajectories/events) */
  TOutputData* OutputData() const { return _fpOutputData; } /**< Output of the simulation */
  long HepRdmSeed() const { return _fHepRdmSeed; } /**< Seed for random numbers, that are managed by 
						    the CLHEP random number class */

 private:
  unsigned long _fN;
  double _fEmin ;
  double _fTmax ;
  string _fRecordMode ;
  TOutputData *_fpOutputData ;
  string _fOutputFileName ;
  long _fHepRdmSeed ;
};

#endif
