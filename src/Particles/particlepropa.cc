/**
   @file    particlepropa.cc
   @author  Eric Armengaud, armengau@in2p3.fr
   @brief   Implementation of the TParticlePropa class. See the .h file
*/

#include "particlepropa.h"
#include <interactiondata.h>
#include <limits>
#include <ctime> 


extern double interactt, deflectt, photodt, pionprodt,  pairpt, decayt;

QUEUE<TParticle*>* TParticlePropa::Propagate(TBasicParam* aBasic) {
  throw  TCrpErr("Propagate of TParticlePropa called!" );
  clock_t ticks;

  if (aBasic->RecordMode() == "Full Trajectories" && _fpUniverse->Type() == UNIVERSE_LARGESCALE)
    aBasic->OutputData()->Add3DTraj(-1,-1,_fInitPosition,_fInitMomentum,-1, -1) ;
  //Allows to separate individual trajectories as well as to give the true initial position
  
  do {
    if (aBasic->RecordMode() == "Full Trajectories") this->Write( aBasic) ;
    if (aBasic->RecordMode() == "Events") this->CheckDetection( aBasic) ;
    
    //_fTimeStep=1.e20;
    //_fNextTimeStep=2.e20;
    
#ifdef DEBUG_OUTPUT
    std::cout<<"\nnucleusId="<<_fChargeNumber*1000 + _fMassNumber<<std::endl;
#endif
#ifdef DEBUG_OUTPUT
    std::cout<<"\nlog10(gamma)="<<log10(_fEnergy/(_fMass*c_squared))<<std::endl;
#endif
 
    

    //TODO check if correct! adaptive stepsize which takes photodisintegration into account
    //TODO This lower bound for TALYS usage should be a parameter
    if(_fMassNumber>1){
      TInteractionData* lpInt = _fpUniverse->InteractionData();
    
      //     double PhotoDisTimeStep = lpInt->GetTimeStep(_fChargeNumber*1000 + _fMassNumber,
      // 						 _fEnergy*(1.+_fRedshift)/1000.,  //GeV
      // 						 _fMass*c_squared/1000.,             //GeV
      // 						 _fRedshift);
   

    
#ifdef DEBUG_OUTPUT
      std::cout<<"AS_TEST: in particle propa before Interact \t_fTimeStep="<<_fTimeStep
	       <<"\t_fNextTimeStep="<<_fNextTimeStep<<std::endl;
#endif
      ticks=clock();
      this->Interact() ; // Change energy, momentum norm
      ticks-=clock();
      interactt+=ticks;
#ifdef DEBUG_OUTPUT
      std::cout<<"AS_TEST: in particle propa after Interact \t_fTimeStep="<<_fTimeStep
	       <<"\t_fNextTimeStep="<<_fNextTimeStep<<"\n\t'*'*'*'*'*'*'*'*'*'*'*'*"<<std::endl;
#endif
 
#ifdef DEBUG_OUTPUT
      //std::cout<<"AS_TEST: PhotoDisTimeStep"<<PhotoDisTimeStep<<std::endl;
#endif
      //_fNextTimeStep=min(_fNextTimeStep,PhotoDisTimeStep);
    }

#ifdef DEBUG_OUTPUT
    std::cout<<"AS_TEST: in particle propa before Deflec \t_fTimeStep="<<_fTimeStep
	     <<"\t_fNextTimeStep="<<_fNextTimeStep<<std::endl;
#endif
    
    ticks=clock();
    this->Deflec() ;   // Change time(steps), position, momentum direction
    ticks-=clock();
    deflectt+=ticks;
#ifdef DEBUG_OUTPUT
    std::cout<<"AS_TEST: in particle propa after Deflec \t_fTimeStep="<<_fTimeStep
	     <<"\t_fNextTimeStep="<<_fNextTimeStep<<std::endl;
#endif
    

  } while (this->CheckEndPropa(aBasic) == 0);
  return &fParticleQueue;
}

bool TParticlePropa::CheckEndPropa(TBasicParam* aBasic) {
  
  bool returnValue = 0 ;

  if (_fEnergy < aBasic->Emin() || _fTime > aBasic->Tmax() ) {
    returnValue=1;
  }
  
  if ( _fpUniverse->Type() == UNIVERSE_ENV1D && _fPosition.x() < _fpUniverse->IntegratorMinTimeStep() * c_light) {
    this->CheckDetection(aBasic); // To avoid final step just to write out the event. 
    returnValue = 1 ;
  }

  return returnValue;
}


void TParticlePropa::CheckDetection( TBasicParam* aBasic) {
  
  if ( _fpUniverse->Observers()->Type() == OBSERVER_POINT ) {
    // 1D test : distance from 0 smaller than next step
    // if (_fPosition.x() <= _fNextTimeStep*c_light ) {
    // Apply correction to reach exactly the origin
    //_fTimeStep = max(_fPosition.x()* inv_c_light, 0.) ;
    //_fPosition += _fTimeStep*c_squared*_fMomentum/_fEnergy ;
    //_fTime += _fTimeStep ;
    // if (_fTimeStep) this->Interact() ;
    _fDistObs=_fPosition.x() ;
    if(_fPosition.x() <= _fpUniverse->IntegratorMinTimeStep() * c_light  ) {
      this->Detect(aBasic) ;
    }
  
    
  } else if ( _fpUniverse->Observers()->Type() == OBSERVER_SMALLSPHERE ) {

        TObservers* lObs = _fpUniverse->Observers() ;
        _fDistObs = DISTOBS_MAX_GPC * Gpc ;
        double lDistObs;
       
        for (int j=0; j<lObs->Nb(); j++) {
	          
	  static  TVector3D lTempPos;
	  static const  double lDeltaX = _fpUniverse->Xmax()-_fpUniverse->Xmin() ;
	  static const  double lDeltaY = _fpUniverse->Ymax()-_fpUniverse->Ymin() ;
	  static const  double lDeltaZ = _fpUniverse->Zmax()-_fpUniverse->Zmin() ;
	  lDistObs= DISTOBS_MAX_GPC * Gpc;

	  for(int ix=-1; ix<=1; ix++){ 
	    for(int iy=-1; iy<=1; iy++){
	      for(int iz=-1; iz<=1; iz++){
		lTempPos.setX( lObs->Positions(j).x()+ lDeltaX * ix) ;
		lTempPos.setY( lObs->Positions(j).y()+ lDeltaY * iy) ;
		lTempPos.setZ( lObs->Positions(j).z()+ lDeltaZ * iz) ;
		lDistObs =min(lDistObs,  (_fPosition - lTempPos).mag() - lObs->Radius()) ;
	      }
	    }
	  }


          if (lDistObs > 0) { // particle is not inside the sphere
    	_fDistObs = min(_fDistObs,lDistObs) ;
	_fDistObs = max(_fDistObs, _fpUniverse->IntegratorMinTimeStep() * inv_c_light) ;

    	if (lDistObs <= lObs->Radius()*OBS_SECURITY_RATIO) { 
    	  // particle position is nearby enough from observer to check if it hits.
    	  double lXP, lDiscrim ;
    	  TVector3D lObsPart = _fPosition - lObs->Positions(j) ;
    	  lXP = lObsPart.dot(_fMomentum) / _fMomentum.mag() ;
    	  lDiscrim = pow(lXP,2) - lObsPart.mag2() + pow(lObs->Radius(),2) ;
    	  if ( lDiscrim >= 0 && lXP <= 0) { 
    	    // Direction defined by momentum hits sphere: particle hits the detector. 
    	    // We correct slightly the position and record
    	    double lTcorr = (-lXP - sqrt(lDiscrim))*(1+OBS_TCORR_INC)* inv_c_light ;
	    static vector<double> lY;
	    lY = this->Y() ;
    	    vector<double> ldYdt; 
	    ldYdt = this->Derivs(lY) ;
    	    for (int i=0; i<6; i++) lY[i] += lTcorr*ldYdt[i] ;
    	    this->YtoXP(lY) ;
    	    _fTime += lTcorr ;
    	    this->Detect(aBasic) ;
    	    _fDistObs = DISTOBS_MAX_GPC * Gpc ;
    	  }  //if true : Direction defined by momentum hits sphere.
    	}  //if true : Particle is close enough to observer to check it it hits.
          }  //if true : Particle not inside sphere
          else{
	  
          }
        }  //loop to identify the distance to the closest observer

        // Bound on next time step due to distance to observers
        _fNextTimeStep = min(_fNextTimeStep,_fDistObs* inv_c_light) ;

  } else if ( _fpUniverse->Observers()->Type() == OBSERVER_LARGESPHERE ) {
 
    // The center of large spheres is taken as initial position.
    // This way we take into account the displacement of the simulation box
    // In this regime the initial position = source = observer obligatorily.
    TObservers* lObs = _fpUniverse->Observers() ;
    _fDistObs = DISTOBS_MAX_GPC * Gpc ;
    double lDistObs ;
    bool lFlag = 0 ;
    for (int j=0; j<lObs->Nb(); j++) {
      lDistObs = (_fPosition-_fInitPosition).mag() - lObs->Radii(j) ;
      _fDistObs = min(_fDistObs, fabs(lDistObs)) ;
      _fDistObs = max(_fDistObs, _fpUniverse->IntegratorMinTimeStep()/ inv_c_light) ;
      if (lDistObs >= 0) lFlag = 0; // part out of the sphere
      else lFlag = 1 ; // particle inside sphere
      if (lFlag != _fIsInsideSpheres[j]) {
	_fIsInsideSpheres[j] = lFlag ;
	this->Detect( aBasic) ;
      }
    }
  } else if (_fpUniverse->Observers()->Type() == OBSERVER_NO) {
    // Nothing !
  } else throw TCrpErr("Error : CheckDetection not implemented in this configuration" );
  
  
}

void TParticlePropa::Detect( TBasicParam* aBasic) const {
  TOutputData* lOutput = aBasic->OutputData() ;
  
  int lType ;
  if (this->Type() == PARTICLE_NUCLEON && _fCharge == 0) lType = 2112 ; // n0
  else if (this->Type() == PARTICLE_NUCLEON && _fCharge == 1) lType = 2212 ; // p+
  else if (this->Type() == PARTICLE_NUCLEUS && _fChargeNumber >= 0 && _fMassNumber>=0) //nucleus
    lType = _fChargeNumber * 1000 + _fMassNumber ;   //Changed by NN 
  else throw TCrpErr("Unknown particle type in Detect.") ;
  
  if (_fpUniverse->Type() == UNIVERSE_ENV1D) {
    lOutput->Add1DEvent(lType,_fInitPosition.mag(),_fInitRedshift,
			_fInitMomentum.mag()*c_light, _fTime,_fEnergy,  _fInitType) ;
  } else if (_fpUniverse->Type() == UNIVERSE_LARGESCALE) {
    lOutput->Add3DEvent(lType,_fInitPosition,_fInitMomentum,_fTime,
			_fPosition,_fMomentum,  _fInitType) ;
  } else 
    throw TCrpErr( "Detection not implemented in this environment" );
  
}

void TParticlePropa::Write( TBasicParam* aBasic) const {
  TOutputData* lOutput = aBasic->OutputData() ;
    
  int lType ;
  if (this->Type() == PARTICLE_NUCLEON && _fCharge == 0) lType = 2112 ; // neutron
  else if (this->Type() == PARTICLE_NUCLEON && _fCharge == 1) lType = 2212 ; // proton
  else if (this->Type() == PARTICLE_NUCLEUS && _fChargeNumber >= 0 && _fMassNumber>=0) //nucleus
    //    lType = _fMassNumber * 1000. + _fChargeNumber ;                     //Changed by NN 
    lType = _fMassNumber + _fChargeNumber * 1000.;                     //Changed by NN 
  else throw TCrpErr("!!!!Unknown particle type in Write.") ;
  
  if (_fpUniverse->Type() == UNIVERSE_ENV1D) {
    lOutput->Add1DTraj(lType,_fTime,_fPosition.mag(),_fEnergy, _fInitType) ;
  } else if (_fpUniverse->Type() == UNIVERSE_LARGESCALE) {
    lOutput->Add3DTraj(lType,_fTime,_fPosition - _fInitPosition,_fMomentum,_fEnergy,  _fInitType) ;
       // lOutput->Add3DTraj(lType,_fTime,_fPosition,_fEnergy,  _fInitType) ;
    // enregistrement du "deplacement", devrait prendre en compte les cbp.
  } else throw TCrpErr("Detection not implemented in this environment");
  
}


void TParticlePropa::Deflec() {
  double  lThisStepSize;
    if(_fNextTimeStep < numeric_limits<double>::min() ) {
      std::cerr<< " WARNING : Propagate : _fNextTimeStep < numeric_limits<double>::min() " << std::endl; 
      std::cerr<< " _fNextTimeStep = MinTimeStep "  << std::endl; 
	_fNextTimeStep=  _fpUniverse->IntegratorMinTimeStep(); 
    }
    lThisStepSize=_fNextTimeStep; 

  if (_fpUniverse->Bflag() && _fCharge != 0) {
    // 3D deflection algorithm
    
    double eps = _fpUniverse->IntegratorEps() ;
    double hdid, hnext, lpropagated;
    double hmax = _fNextTimeStep;// min(_fpUniverse->InteractionData()->InteractionTimeStep(),_fDistObs) ; // Sollte _fNextTimeStep sein?? 
    lpropagated=0; 
    // Upper value on the time step: 1) interaction step, 
    //   2) distance to nearest observer
    if (_fpUniverse->Integrator() == "Cash-Karp RK") {
      static int deFlag=0;
      if(deFlag == 0){ 
	deFlag=1;
	std::cout<< "In deflec()" << std::endl;
      }

      static     vector<double> lY;
      static     vector<double> ldYdt;
      static    vector<double> lYscal(6,0) ;
      lY = this->Y() ;
      ldYdt = this->Derivs(lY ) ;
      //Old lYscal definition. 

      for (int i=0; i<6; i++) lYscal[i] = fabs(lY[i]) + fabs(_fNextTimeStep*(ldYdt)[i]) ;
      while(_fNextTimeStep > _fpUniverse->IntegratorMinTimeStep() ){// This while grants that the CR is propagated over lpropagated.
	

	//Proposed new definition of lYscal definition. 
	//for (int i=0; i<6; i++) lYscal[i] = abs(_fNextTimeStep*(ldYdt)[i]) ;
	
	lY = this->Y() ;
        ldYdt= this->Derivs(lY ) ;
	this->RKqs(lY, ldYdt, _fTime, min(_fNextTimeStep, 1 *(Mpc * inv_c_light) ), eps, lYscal, hdid, hnext, hmax) ;
	this->YtoXP(lY) ;
	//this->RKqs(lY, ldYdt, _fTime, _fNextTimeStep, eps, lYscal, hdid, hnext, hmax) ;
	lpropagated += hdid;
	_fNextTimeStep=lThisStepSize - lpropagated;

      }

      /* BS integrator not available yet!
	 } else if (_fpUniverse->Integrator() == "Bulirsch-Stoer") {
	 vector<double> lY = this->Y() ;
	 vector<double> ldYdt = this->Derivs(lY ) ;
	 vector<double> lYscal(6,0) ;
	 for (int i=0; i<6; i++) lYscal[i] = abs(lY[i]) + abs(_fNextTimeStep*(ldYdt)[i]) ;
	 this->BulirschStoer( lY, ldYdt, _fTime, _fNextTimeStep, eps, 
	 lYscal, hdid, hnext, hmax) ;
	 this->YtoXP(lY) ;
      */
    } else {
      throw TCrpErr( "Integrator not implemented. Exiting" );
    }

    _fTime += lpropagated ;
  }
  
  // Straight line propagation
  _fTimeStep = lThisStepSize ;
  _fPosition += _fNextTimeStep*c_squared*_fMomentum/_fEnergy ; //linear step to reach exactly _fNextTimeStep
  _fTime += _fNextTimeStep;


  if (_fpUniverse->Type() == UNIVERSE_LARGESCALE) {
    // Periodic boudary conditions applied :
    // Position remains in the box, but initial position is increased
    while (_fPosition.x() < _fpUniverse->Xmin()) {
      double lDelta = _fpUniverse->Xmax()-_fpUniverse->Xmin() ;
      _fPosition.setX(_fPosition.x()+lDelta) ;
      _fInitPosition.setX(_fInitPosition.x()+lDelta) ;
    }
    while (_fPosition.x() >= _fpUniverse->Xmax()) {
      double lDelta = _fpUniverse->Xmax()-_fpUniverse->Xmin() ;
      _fPosition.setX(_fPosition.x()-lDelta) ;
      _fInitPosition.setX(_fInitPosition.x()-lDelta) ;
    }
    while (_fPosition.y() < _fpUniverse->Ymin()) {
      double lDelta = _fpUniverse->Ymax()-_fpUniverse->Ymin() ;
      _fPosition.setY(_fPosition.y()+lDelta) ;
      _fInitPosition.setY(_fInitPosition.y()+lDelta) ;
    }
    while (_fPosition.y() >= _fpUniverse->Ymax()) {
      double lDelta = _fpUniverse->Ymax()-_fpUniverse->Ymin() ;
      _fPosition.setY(_fPosition.y()-lDelta) ;
      _fInitPosition.setY(_fInitPosition.y()-lDelta) ;
    }

    while (_fPosition.z() >= _fpUniverse->Zmax()) {
      double lDelta = _fpUniverse->Zmax()-_fpUniverse->Zmin() ;
      _fPosition.setZ(_fPosition.z()-lDelta) ;
      _fInitPosition.setZ(_fInitPosition.z()-lDelta) ;
    }

    while (_fPosition.z() < _fpUniverse->Zmin()) {
      double lDelta = _fpUniverse->Zmax()-_fpUniverse->Zmin() ;
      _fPosition.setZ(_fPosition.z()+lDelta) ;
      _fInitPosition.setZ(_fInitPosition.z()+lDelta) ;
    }
  }
}

inline TVector3D TParticlePropa::LorentzForce(TVector3D aPos, TVector3D aSpeed) const {
  static  TVector3D lField;
  lField = _fpUniverse->MagField()->getField(aPos) ;
  static  TVector3D lForce;
  lForce = _fCharge*aSpeed.cross(lField);
  return lForce;
}

vector<double> TParticlePropa::Y() const {
 static vector<double> Y(6,0) ;
  Y[0] = _fPosition.x() ;
  Y[1] = _fPosition.y() ;
  Y[2] = _fPosition.z() ;
  Y[3] = _fMomentum.x() ;
  Y[4] = _fMomentum.y() ;  
  Y[5] = _fMomentum.z() ;
  return Y;
}

void TParticlePropa::YtoXP(vector<double> aY) {
  
  _fPosition.set(aY[0],aY[1],aY[2]) ;
  _fMomentum.set(aY[3],aY[4],aY[5]) ;
  
  if (KILL_NUM_DISS == 1) {
    // rustine pour tuer les dissipations numeriques d'energie...
    _fMomentum.setMag(_fEnergy* inv_c_light) ;
  }

}

vector<double> TParticlePropa::Derivs(vector<double> Y) const {
  // 6-D derivative. array sizes must be 6 only here (not specified as args).
  // Use of standard Hamiltonian formalism
  // Y=(pos,momentum) Derivs=(speed,force)
static  vector<double> dY(6,0) ;
 static TVector3D lSpeed, lPos, lForce;
  double normp = sqrt(pow(Y[3],2)+pow(Y[4],2)+pow(Y[5],2)) ;
  dY[0] = c_light*Y[3]/normp ;
  dY[1] = c_light*Y[4]/normp ;
  dY[2] = c_light*Y[5]/normp ;
  //  TVector3D lSpeed(dY[0],dY[1],dY[2]) ;
  lSpeed[0]=dY[0];
  lSpeed[1]=dY[1];
  lSpeed[2]=dY[2];
  //  TVector3D lPos(Y[0],Y[1],Y[2]) ;
  lPos[0]=Y[0];
  lPos[1]=Y[1];
  lPos[2]=Y[2];
  //  TVector3D lForce = LorentzForce(lPos, lSpeed ) ;
 lForce = LorentzForce(lPos, lSpeed ) ;
  dY[3] = lForce.x() ;
  dY[4] = lForce.y() ; 
  dY[5] = lForce.z() ;
  return dY;
}

void TParticlePropa::RKqs( vector<double>& Y, vector<double> dYdt, double t, double htry, double eps, vector<double> Yscal, double& hdid, double& hnext, double hmax) {
	
  double errmax, htemp, tnew ;
  vector<double> Yerr(6,0) ;
  vector<double> Ytemp(6,0) ;
  double h = htry ;
  // Changement : def. de hmin inclut DistObs!
  double hmin = min(_fpUniverse->IntegratorMinTimeStep(),_fDistObs/c_light) ;

  while(1) {
    // Loop to select the best stepsize as possible, within the simulation constraints.
    this->RKCK(Y,dYdt,h,Ytemp,Yerr) ;
    errmax = 0.0 ;
    for (int i=0; i<6; i++)
      // Changement : essai de meilleure contrainte sur Yscal/Yerr!
      if (Yscal[i] > numeric_limits<double>::min() ) errmax=max(errmax,fabs(Yerr[i]/Yscal[i])) ; 
    /*
      1. Yscal estimated by Euler integration Yscal=Y0 + dxdy * timestep
      2. Yerr estimated fith order error of the RKCK integrator. 
      3. Yscal proportional to estimated stepsize for this step (_fnextTimeStep). This should change with h. 
    */

    errmax /= eps ;
    // Changement : un seul break, avant de modifier h!
    if (errmax <= 1 || h <= 1.001*hmin) break ;
    htemp = RK_SAFETY*h*pow(errmax,RK_PSHRNK) ;
    h=max(htemp,0.1*h); // h > 0
    // Changement : contrainte sur h en plus!
    h = max(hmin,min(h,hmax)) ;
    // if (h <= hmin) break ; // New w.r.t. Num. recipes : due to exterior constraints on the stepsize.
    tnew=t+h ;
    if (tnew == t) throw TCrpErr("stepsize underflow in rkqs" );
  }
  
  if (errmax > RK_ERRCON) {
    hnext = RK_SAFETY*h*pow(errmax,RK_PGROW) ;
  } else {
    hnext = 5*h ;
  }
  hnext = max(hmin,min(hnext, hmax)) ;
  hdid = h ;
  for (int i=0; i<6; i++) Y[i] = Ytemp[i] ;
  
}

/*
 
  void TParticlePropa::BulirschStoer( vector<double>& Y, vector<double> dYdt, double t, double htry, double eps, vector<double> Yscal, double& hdid, double& hnext, double hmax) {
	 
  //routine derived from f77 code
  // NOT WORKING YET !!!!!!
	 
  vector<double> lYerr(6,0) ;
  vector<double> lY0 = Y;
  vector<double> ldY0 = dYdt;
  vector<double> lYtemp(6,0), lYseq(6,0) ;
  double t0 = t ;
  double h,xest,a0,errmax ;
  int i = 0 ;
  int i_err = 0 ;
  static int Nseq[BS_IMAX] = {2,4,6,8,12,16,24,32,48,64,96} ;
	 
  double hmin = _fpUniverse->IntegratorMinTimeStep() ;
	 
	 
  h = htry ;
  while (i_err == 0 && i < BS_IMAX) {
  lYseq = this->Mmid(lY0, ldY0, t0, h, Nseq[i], ) ;
  a0 = h/Nseq[i] ;
  xest = a0*a0 ;
  this->Rzextr(i,xest,lYseq,Y,lYerr) ;
  errmax = 0 ;
  for (int j=0; j<6; j++) errmax = max(errmax,abs(lYerr[j]/Yscal[j])) ;
  errmax /= eps ;
  if (errmax < 1) {
  t += h ;
  hdid = h ;
  if (i == BS_NUSE-1) hnext = h*BS_SHRINK ;
  else if ( i == BS_NUSE-2 ) hnext = h*BS_GROW ;
  else hnext = (h*Nseq[BS_NUSE-2])/Nseq[i] ;
  i_err = 1 ;
  }
  i += 1 ;
  if (i >= BS_IMAX && i_err == 0) {
  h /= 10. ;
  if (t+h == t) cerr << "stepsize underflow in BS" << endl;
  i = 1 ;
  }
  }
	 
  }
 
  void TParticlePropa::Rzextr(int iest, double xest, vector<double> Yest, vector<double>& Yz, vector<double>& dY) {
	 
  double yy,ddy,b,c,b1,v ;
  int m1 ;
  //  double D[BS_NMAX][BS_NCOL] ;
  //  double x[BS_IMAX] ;
  double fx[BS_NCOL] ;
	 
  _fBS_x[iest] = xest ;
  if (iest == 0) {
  Yz = Yest ;
  dY = Yest ;
  for (int j=0; j<6; j++) _fBS_D[j][0] = Yest[j] ;
  } else {
  m1 = min(iest,BS_NUSE-1) ;
  for (int k=0; k<m1-1; k++) fx[k+1] = _fBS_x[iest-k-1]/xest ;
  for (int j=0; j<6; j++) {
  yy = Yest[j] ;
  v = _fBS_D[j][0] ;
  c = yy ;
  _fBS_D[j][0] = yy ;
  for (int k=1; k<m1; k++) {
  b1 = fx[k]*v ;
  b = b1-c ;
  ddy = v ;
  if (b != 0) {
  b = (c-v)/b ;
  ddy = c*b ;
  c = b1*b ;
  }
  if (k != m1) v = _fBS_D[j][k] ;
  _fBS_D[j][k] = ddy ;
  yy += ddy ;
  }
  dY[j] = ddy ;
  Yz[j] = yy ;
  }
  }
	 
  }
 
  vector<double> TParticlePropa::Mmid(vector<double> Y, vector<double> dYdt, double xs, double htot, int Nstep, ) {
	 
  vector<double> Ym = Y ;
  vector<double> Yn(6,0) ;
  vector<double> Yout(6,0) ;
  double h,x,swap,h2;
	 
  h = htot/Nstep ;
  for (int i=0; i<6; i++) Yn[i] = Y[i] + h*dYdt[i] ;
  x = xs+h ;
	 
  Yout = this->Derivs(Yn, ) ;
  h2 = 2*h ;
  for (int n=1; n<Nstep; n++) {
  for (int i=0; i<6; i++) {
  swap = Ym[i]+h2*Yout[i] ;
  Ym[i] = Yn[i] ;
  Yn[i] = swap ;
  }
  x += h ;
  Yout = this->Derivs(Yn, ) ;
  }
  for (int i=0; i<6; i++) Yout[i] = 0.5*(Ym[i]+Yn[i]+h*Yout[i]) ;
	 
  return Yout ;
  }
 
*/

void TParticlePropa::RKCK( vector<double> Y, vector<double> dYdt, double h, vector<double>& Yout, vector<double>& Yerr) {
  // 5th order Cash-Karp RK method
  //static double a2=0.2,a3=0.3,a4=0.6,a5=1.0,a6=0.875 ;
  // coefficients a_i are not needed as long as Derivs does not depend on time
  static const double b21=0.2,
    b31=3.0/40.0,b32=9.0/40.0,b41=0.3,b42=-0.9,b43=1.2,
    b51=-11.0/54.0,b52=2.5,b53=-70.0/27.0,b54=35.0/27.0,
    b61=1631.0/55296.0,b62=175.0/512.0,b63=575.0/13824.0,
    b64=44275.0/110592.0,b65=253.0/4096.0,c1=37.0/378.0,
    c3=250.0/621.0,c4=125.0/594.0,c6=512.0/1771.0,
    dc5=-277.0/14336.0 ;
  static const double dc1=c1-2825.0/27648.0,dc3=c3-18575.0/48384.0,
    dc4=c4-13525.0/55296.0,dc6=c6-0.25 ;
  
  //#define USE_STL_ALGO
#ifdef USE_STL_ALGO 
  
  static vector<double> ak2;
  static vector<double> ak3;
  static vector<double> ak4;
  static vector<double> ak5;
  static vector<double> ak6;
  
  static vector<double> Ytemp(6) ;
  static vector<double> Y2(6);
  static vector<double> Y3(6);
  static vector<double> Y4(6);
  
  A_EQ_B_MULT_K ( Y2, dYdt, b21*h ) ;
  A_EQ_B_PLUS_C ( Ytemp, Y, Y2 );
  ak2 = this->Derivs(Ytemp ) ;
  A_EQ_B_MULT_K ( Y2 , ak2 , b32*h ) ;
  A_EQ_B_MULT_K ( Y3 , dYdt, b31*h ) ;
  A_EQ_B_PLUS_C ( Y4 , Y2 , Y3 ) ;
  A_EQ_B_PLUS_C ( Ytemp, Y, Y4 ) ;
  ak3 = this->Derivs(Ytemp ) ;
  A_EQ_B_MULT_K ( Y2 , dYdt , b41*h );
  A_EQ_B_MULT_K ( Y3 , ak2  , b42*h );
  A_EQ_B_PLUS_C ( Y4 , Y2   , Y3 ) ;
  A_EQ_B_MULT_K ( Y2 , ak3  , b43*h ) ;
  A_EQ_B_PLUS_C ( Y3 , Y2   , Y ) ;
  A_EQ_B_PLUS_C ( Ytemp, Y3 , Y4 );
  ak4 = this->Derivs(Ytemp ) ;
  A_EQ_B_MULT_K ( Y2 , dYdt , b51*h );
  A_EQ_B_MULT_K ( Y3 , ak2  , b52*h );
  A_EQ_B_PLUS_C ( Y4 , Y2 , Y3 ) ;
  A_EQ_B_MULT_K ( Y2 , ak3  , b53*h );
  A_EQ_B_PLUS_C ( Y3 , Y4 , Y2 ) ;
  A_EQ_B_MULT_K ( Y2 , ak4  , b54*h );
  A_EQ_B_PLUS_C ( Y4 , Y3 , Y2 ) ;
  A_EQ_B_PLUS_C ( Ytemp , Y4, Y );
  ak5 = this->Derivs(Ytemp ) ;
  
  for (int i=0; i<6; i++) Ytemp[i]=Y[i]+h*(b61*dYdt[i]+b62*ak2[i]+b63*ak3[i]+b64*ak4[i]+b65*ak5[i]) ;
  ak6 = this->Derivs(Ytemp ) ;
  for (int i=0; i<6; i++) Yout[i] = Y[i]+h*(c1*dYdt[i]+c3*ak3[i]+c4*ak4[i]+c6*ak6[i]) ;
  for (int i=0; i<6; i++) Yerr[i] = h*(dc1*dYdt[i]+dc3*ak3[i]+dc4*ak4[i]+dc5*ak5[i]+dc6*ak6[i]) ;
  
#else
  
static  vector<double> ak2(6,0) ;
static  vector<double> ak3(6,0) ;
static  vector<double> ak4(6,0) ;
static  vector<double> ak5(6,0) ;
static  vector<double> ak6(6,0) ;
static  vector<double> Ytemp(6,0) ;
  for (int i=0; i<6; i++) Ytemp[i]=Y[i]+b21*h*dYdt[i] ;
  ak2 = this->Derivs(Ytemp ) ;
  for (int i=0; i<6; i++) Ytemp[i]=Y[i]+h*(b31*dYdt[i]+b32*ak2[i]) ;
  ak3 = this->Derivs(Ytemp ) ;
  for (int i=0; i<6; i++) Ytemp[i]=Y[i]+h*(b41*dYdt[i]+b42*ak2[i]+b43*ak3[i]) ;
  ak4 = this->Derivs(Ytemp ) ;
  for (int i=0; i<6; i++) Ytemp[i]=Y[i]+h*(b51*dYdt[i]+b52*ak2[i]+b53*ak3[i]
					   +b54*ak4[i]) ;
  ak5 = this->Derivs(Ytemp ) ;
  for (int i=0; i<6; i++) Ytemp[i]=Y[i]+h*(b61*dYdt[i]+b62*ak2[i]+b63*ak3[i]
					   +b64*ak4[i]+b65*ak5[i]) ;
  ak6 = this->Derivs(Ytemp ) ;
  for (int i=0; i<6; i++) Yout[i] = Y[i]+h*(c1*dYdt[i]+c3*ak3[i]
					    +c4*ak4[i]+c6*ak6[i]) ;
  for (int i=0; i<6; i++) Yerr[i] = h*(dc1*dYdt[i]+dc3*ak3[i]+dc4*ak4[i]+dc5*ak5[i]+dc6*ak6[i]) ;
  
#endif
  
}

void TParticlePropa::RedshiftLoss() {
  //std::cout<<"Start RedshiftLoss()."<<std::endl;

  static int testflag=0; 
 if(testflag==0){
  testflag=1;
  std::cout<<"\t->RedshiftLoss called" << std::endl; 
  //std::cout<< "Achtung modified!!! " << std::endl; 
 }

  if (_fpUniverse->Type() == UNIVERSE_ENV1D) {
    // Using the relations dE/dz=E/(1+z) and dz/dt=H_0(1+z)E(z)
    double lLoss = _fpUniverse->H0()*_fTimeStep*sqrt(_fpUniverse->OmegaLambda()+_fpUniverse->OmegaM()*pow(1+_fRedshift,3)) ;
    _fEnergy *= (1-lLoss) ;
    _fMomentum.setMag(_fEnergy * inv_c_light) ;
  }

}
