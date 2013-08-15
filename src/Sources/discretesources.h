/**
   @file   discretesources.h
   @author Eric Armengaud, armengau@in2p3.fr
   @brief  Class describing discrete sources (position, spectrum)
*/

#ifndef _DISCRETESOURCES_H_
#define _DISCRETESOURCES_H_

#include <vector>

#include "sources.h"
#include "xmlparam.h"
#include "crp_err.h"


#include <fstream>
#include "CLHEP/Random/RandFlat.h"

using namespace std;

/**
   @class TDiscreteSources
   @brief Class describing discrete sources: a finite number of sources are built once of all at the beginning of the simulation, and one of these sources is choosen at random for each injected particle during the simulation.

   Discrete sources can be either specified from a list of positions in the configuration file, or they can be drawn at random at the beginning of the simulation from a source density grid. In this case, the configuration file must indicate the name of the density file and the number of sources that will be drawn.
*/

class TDiscreteSources : public TSources, public TXmlParam {
 public:
  TDiscreteSources(const char*) ;
  ~TDiscreteSources() {}

  unsigned int Nb() const { return _fNumber; }
  /**< Number of sources */
  TVector3D Positions(int i) const { return _fPositions.at(i); }
  /**< Finite list of source positions */
  double EcutList(int i) const { return _fEcutList.at(i); }
  /**< List of maximum energies (case of a power-law spectrum) or injection energies (case of a monochromatic spectrum) for all the sources. NOT YET IMPLEMENTED. */
  double AlphaList(int i) const { return _fAlphaList.at(i); }
  /**< List of spectral indices for each individual source in the case of a power-law spectrum. These indices are drawn at the beginning of the simulation from a top-hat distribution in the range [Alpha()-SigAlpha() , Alpha()+SigAlpha()] */
  double SigAlpha() const { return _fSigAlpha; }
  /**< Half-width of the distribution of spectral indices for the sources. See AlphaList(). This value is set by default to 0. */

 private:
  unsigned int _fNumber;
  double _fSigAlpha ;
  vector<TVector3D> _fPositions ;
  vector<double> _fAlphaList ;
  vector<double> _fEcutList ;
  vector<double> _fEminList ;

};

#endif
