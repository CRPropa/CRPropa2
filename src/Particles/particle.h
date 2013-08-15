
/**
   @file    particle.h
   @author  Tristan Beau, beau@in2p3.fr
   @brief   Generic class for particles
*/

#ifndef _PARTICLE_H_
#define _PARTICLE_H_

#include "typedclass.h"
#include "universe.h"
#include "basicparam.h"
#include "vector3d.h"
#include "defqueue.h"
#include <iostream>
#include <vector>

using namespace std;

/**
   @class TParticle
   @brief Generic virtual class for any particle

*/

class TParticle : public TTypedClass {
 public:
  virtual ~TParticle() { }

  double Mass() const { 
    std::cout<<" Non overloaded Mass()"<<std::endl;
    return _fMass; }
  /**< Particle mass */
  double Charge() const { return _fCharge; }
  /**< Particle charge */
  TVector3D Position() const { return _fPosition; }
  /**< Position = vector r */
  TVector3D Momentum() const { return _fMomentum; }
  /**< Momentum = vector p */
  double Energy() const { return _fEnergy; }
  /**< Energy = c_light * Momentum */
  double Redshift() const { return _fRedshift; }
  /**< Particle redshift. Set to 0 in 3-dimensionnal simulations. */
  void ComputeRedshift() ; 
  /**< Compute the redshift from the particle position and the cosmological parameters. */

  
  //Changed by NN
  int MassNumber() const { return _fMassNumber; }
  /**< Returns the particle's mass number. */
  int ChargeNumber() const { return _fChargeNumber; }
  /**< Returns the particle's charge number. */
#ifdef EXTENDED_ROOT_DATA_STRUCTURE 
  long unsigned int GetEveId() const { return _fEveId; }
#endif
  
  virtual QUEUE<TParticle*>* Propagate(TBasicParam*)
    { return 0; }
  static QUEUE<TParticle*> fParticleQueue; 
  /**< List of pointers to the secondary particles which still have to be treated before we generate a new particle from the sources. When a secondary is generated during propagation, it is added to that Queue. */

 protected:
#ifdef EXTENDED_ROOT_DATA_STRUCTURE 
  unsigned long int _fEveId;  
  unsigned long int _fMotherId;
  unsigned long int _fDaugtherId;
#endif
  
  double _fMass;
  double _fCharge;
  int _fMassNumber;      //changed by NN
  int _fChargeNumber;    //changed by NN
  double _fEnergy ;
  TVector3D _fPosition;
  TVector3D _fMomentum;
  double _fRedshift;
  double _fInitRedshift;
  TUniverse* _fpUniverse;
    std::vector<double> da;

};

#endif
