/**
   @file    particlepropa.h
   @author  Eric Armengaud, armengau@in2p3.fr
   @brief   Class of particles with propagation properties
*/

#ifndef _PARTICLEPROPA_H_
#define _PARTICLEPROPA_H_

#include <fstream>
#include <math.h>
#include <vector>
#include <functional>
#include <cfloat>
#include <limits.h>
#include "particle.h"

#include "universe.h"
#include "vector3d.h"
#include "basicparam.h"
#include "env1d.h"
#include "units.h"
#include "algo_vec.h"



#define RK_SAFETY 0.9
/**< Parameter for the 5th order Runge-Kutta integrator */
#define RK_PGROW -0.2
/**< Parameter for the 5th order Runge-Kutta integrator */
#define RK_PSHRNK -0.25
/**< Parameter for the 5th order Runge-Kutta integrator */
#define RK_ERRCON 1.89e-4
/**< Parameter for the 5th order Runge-Kutta integrator */
/* BS integrator parameters
  #define BS_NMAX 10
  #define BS_IMAX 11
  #define BS_NUSE 7
  #define BS_ONE 1
  #define BS_SHRINK 0.95
  #define BS_GROW 1.2
  #define BS_NCOL 7
*/

#define DISTOBS_MAX_GPC 1.
/**< Initial value for the the distance to observers. Once a particle has been detected by a small sphere observer, its parameter DistObs, used to test detection, is reset to DISTOBS_MAX_GPC. */

#define KILL_NUM_DISS 1
/**< If set, remove numerical dissipation effects on energy in integrator */

using namespace std;

/**
   @class TParticlePropa
   @brief Virtual class describing particles whose propagation is followed "step by step", ie. charged particles.

*/

class TParticlePropa : public TParticle {

 public:
  double Time() const { return _fTime; }
  /**< Propagation time, set to 0 at the beginning. */
  double TimeStep() const { return _fTimeStep; }
  /**< Effective timestep at a given step */
  double NextTimeStep() const { return _fNextTimeStep; } 
  /**< In 3D, it represents the initial timestep that will be used in the integrator during the next step. */
  double DistObs() const { return _fDistObs; }
  /**< Distance to the nearest small sphere observer. The time steps are computed so that they never run beyond this value. */
  TVector3D InitPosition() const { return _fInitPosition; }
  /**< Initial particle position. In 3D, because of periodic boundary conditions, the initial position is shifted each time the particle changes of "box". As a consequence, the particle position is always in the range of the simulation box limits, but its initial position can go well beyond. */
  TVector3D InitMomentum() const { return _fInitMomentum; }
  /**< Initial particle momentum. Contains initial energy as well as injection direction. */
  bool IsInsideSpheres(int i) const { return _fIsInsideSpheres.at(i); }
  /**< Flags to know wether the particle is inside each of the large spherical observers. A particle is recorded if this flag changes during a given step. */

  virtual  QUEUE<TParticle*>* Propagate(TBasicParam*) ;
  /**< Driver routine. Main CRPropa routine, implementing the loop on propagation steps. */
  void Deflec() ;
  /**< At each step, changes position and momentum direction. This is trivially achieved in 1 dimension, and requires the integration of the Lorentz force in 3 dimensions. */
  virtual void Interact() = 0 ;
  /**< At each step, change energy and momentum norm. This takes into account all interactions, and generates secondaries if necessary. */
  bool CheckEndPropa(TBasicParam*) ;
  /**< Determines when the propagation finishes. This is the case when the energy goes below Emin, or when the particle has been propagated during a time larger than Tmax. */
  void RedshiftLoss() ;
  /**< In 1 dimension, changes the energy and momentum due to redshift effects at each step. Uses the cosmological parameters of the simulation. */

  void Write(TBasicParam*) const ;
  /**< To write out a full trajectory. */
  void CheckDetection(TBasicParam*) ;
  /**< To check if there is detection by an observer, and record an "event" if necessary. */
  void Detect(TBasicParam*) const ;   
  /**< To print an "event" into the output if there is detection. */

  TVector3D LorentzForce(TVector3D, TVector3D) const ; 
  /**< Compute Lorentz force from position and speed vectors */
  vector<double> Y() const ;
  /**< In 6D phase space, Y = (position,momentum) */
  void YtoXP(vector<double>) ;
  /**< Phase space <-> real space conversion. */
  vector<double> Derivs(vector<double>) const ; 
  /**< Compute the time derivative dY/dt in phase space, using the Lorentz force. */
  void RKCK(vector<double>, vector<double>, 
	    double, vector<double>&, vector<double>&) ;
  /**< 5th order RK integrator. From the Num. Recipes. */
  void RKqs(vector<double>&, vector<double>, double, double, 
	    double, vector<double>, double&, double&, double) ; 
  /**< Adaptative stepsize driver for RK integrator. From the Num. Recipes */
virtual void PairProduction( double lTimeStep)=0 ;
  /* BS integrator is not working yet !
     void BulirschStoer(const TUniverse*, vector<double>&, vector<double>, double, double, double, vector<double>, double&, double&, double) ;
     void Rzextr(int, double, vector<double>, vector<double>&, vector<double>&) ;
     vector<double> Mmid(vector<double>, vector<double>, double, double, int, const TUniverse*) ;
  */



 protected:
  double _fTime;
  double _fTimeStep;
  double _fNextTimeStep;
  double _fDistObs;
  bool _sDetected;
  vector<bool> _fIsInsideSpheres ;
  vector<bool> _fWasInsideSpheres ;
  TVector3D _fInitPosition;
  TVector3D _fInitMomentum;

  double _fInitRedshift;

  int _fInitType;
  /* BS integrator is not working yet!
     // specific to BS integrator (rational extrapolation) - 
     // better than an "extern"...
     double _fBS_D[BS_NMAX][BS_NCOL] ;
     double _fBS_x[BS_IMAX] ;
  */

};

#endif
