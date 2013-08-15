/**
   @file    basicpinteractions.h
   @author  Eric Armengaud, armengau@in2p3.fr
   @brief   Class describing interaction data as modelled in the F77 code
*/

#ifndef _BASICPINTERACTIONS_H_
#define _BASICPINTERACTIONS_H_

#include "interactiondata.h"
#include "xmlparam.h"
#include "crp_err.h"

#include <fstream>
#include <math.h>
#include "units.h"

// TODO : METTRE DES WARNINGS QD EMIN/MAX DEPASSENT LE DOMAINE TABULE...

#define LOW_ENERGY_BASICPINTERACTIONS 9.*EeV
/**< Used when a particle's energy goes below the tabulated range, which is 10 EeV in this case. */

using namespace std;

/**
   @class TBasicPInteractions
   @brief Stand-alone interaction tables for the proton and neutron pion production on the CMB and for the proton pair production on the CMB.

   This is the version of interactions inherited from the former Fortran code. Redshift effects, pair and pion productions can be switched on/off by the configuration file. No secondary generation is possible. We recommand to use the more complete TSophiaInteractions class.

*/

class TBasicPInteractions : public TInteractionData, public TXmlParam {
 public:
  TBasicPInteractions(const char*) ;
  ~TBasicPInteractions() ;

  double L_pair(int i) const { return _fL_pair[i]; }
  /**< Energy loss length tabulated for pair production on the CMB. */ 
  double E_part(int i) const { return _fE_part[i]; }
  /**< Energy grid used for all tabulations. */
  double dEtabbin() const { return _fdEtabbin; }
  /**< Logarithmic bin of the energy grid. */
  double inc_p(int i) const { return _finc_p[i]; }
  double inc_n(int i) const { return _finc_n[i]; }
  double Eminp(int i) const { return _fEminp[i]; }
  double Eminn(int i) const { return _fEminn[i]; }
  double lossp(int i, int j) const { return _flossp[i][j]; }
  double lossn(int i, int j) const { return _flossn[i][j]; }
  double Lpp(int i, int j) const { return _fLpp[i][j]; }
  double Lnn(int i, int j) const { return _fLnn[i][j]; }
  double lossp_tot(int i) const { return _flossp_tot[i]; }
  /**< Total interaction rate for protons by pair production on the CMB. */
  double lossn_tot(int i) const { return _flossn_tot[i]; }
  /**< Total interaction rate for neutrons by pair production on the CMB. */

  bool PairProdFlag() const { return _fPairProdFlag; }
  /**< If set, pair production is taken into account. Set to 1 by default. */
  bool PionProdFlag() const { return _fPionProdFlag; }
  /**< If set, pion production is taken into account. Set to 1 by default. */
  bool RedshiftFlag() const { return _fRedshiftFlag; }
  /**< If set, redshift losses are taken into account. Set to 1 by default. */

  RandGeneral* RandDistriNeutronDecay() const { return _fpRandDistriNeutronDecay; }
  /**< Energy distribution of the outcoming electron when a neutron decay takes place. */

 protected:
  double _fL_pair[81] ; // parameters in the common sections
  double _fE_part[81] ; // pairloss, Egrid, pionI,II,III
  double _fdEtabbin ;   // of mgrid_far.f
  double _finc_p[81],_finc_n[81] ;
  double _fEminp[81],_fEminn[81] ;
  double _flossp[81][81],_flossn[81][81] ;
  double _fLpp[81][81],_fLnn[81][81] ;
  double _flossp_tot[81],_flossn_tot[81] ;

  bool _fPairProdFlag ;
  bool _fPionProdFlag ;
  bool _fRedshiftFlag ;

  RandGeneral* _fpRandDistriNeutronDecay ;

};

#endif
