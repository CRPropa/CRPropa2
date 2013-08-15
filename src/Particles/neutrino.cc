

/**
   @file    neutrino.cc
   @author  Eric Armengaud, armengau@in2p3.fr
   @brief   Implementation of the TNeutrino class. See the .h file
*/

#include "neutrino.h"
#include <exception>

TNeutrino::TNeutrino(TUniverse* aUniv)  {
  SetType(PARTICLE_NEUTRINO);
  _fpUniverse = aUniv;
  throw TCrpErr("neutrino source not implemented") ;
  
}

TNeutrino::TNeutrino(int aFlavor, TVector3D aPosition, TVector3D aMomentum, double aSourceEnergy, TUniverse* aUniv, int aInitType, double aInitTime) {
  SetType(PARTICLE_NEUTRINO);
  _fFlavor = aFlavor ;
  _fCharge = 0 ;
  _fMass = 0 ;
  _fInitPosition = aPosition ;
  _fPosition = aPosition ;
  _fMomentum = aMomentum ;
  _fEnergy = _fMomentum.mag() * c_light ;
  _fInitialEnergy = _fEnergy ;
  _fInitType=aInitType;
  _fSourcePosition = TVector3D() ;
  _fSourceEnergy = aSourceEnergy ;
  _fInitType=aInitType;
  _fpUniverse = aUniv;
  _fInitTime = aInitTime;
    da = _fpUniverse->DistanceArray();
}

TNeutrino::TNeutrino(int aFlavor, TVector3D aPosition, TVector3D aMomentum, TVector3D aSourcePosition, double aSourceEnergy, TUniverse* aUniv, int aInitType, double aInitTime) {
  SetType(PARTICLE_NEUTRINO);
  _fFlavor = aFlavor ;
  _fCharge = 0 ;
  _fMass = 0 ;
  _fSourcePosition = aSourcePosition ;
  _fInitPosition = aPosition ;
  _fPosition = aPosition ;
  _fMomentum = aMomentum ;
  _fEnergy = _fMomentum.mag() * c_light ;
  _fInitialEnergy = _fEnergy ;
  _fSourceEnergy = aSourceEnergy ;
  _fInitType=aInitType;
  _fpUniverse = aUniv;
  _fInitTime = aInitTime;
    da = _fpUniverse->DistanceArray();
}

QUEUE<TParticle*>* TNeutrino::Propagate(TBasicParam* aBasic) {
  

#ifdef DEBUG_OUTPUT
  std::cout<< "Start of neutrino propagation " << std::endl;
#endif
  //#define NU_TEST
#ifdef NU_TEST
  _fInitPosition.set(3 * Mpc, 3 * Mpc, 3 * Mpc);
  _fPosition=_fInitPosition;
  _fMomentum.setRThetaPhi( _fEnergy* inv_c_light, 3.141592653 * 0.1475846, 3.141592653 * 3/2. ) ;
#endif
  vector<double> lCheck = this->CheckDetection(aBasic) ;
#ifdef NU_TEST
  std::cout<< "Neutrino Propagation modified!!!" << std::endl;
  std::cout<< "lCheck.size() " << lCheck.size() << std::endl;
  for(int i=0; i< lCheck.size() ; i++) std::cout<< "lCheck[" << i << "] " << lCheck[i] << std::endl;
#endif
  // CheckDetection also corrects for the neutrino position if there has been a detection
  if (lCheck.size() > 0) {
    if (_fpUniverse->Type() == UNIVERSE_ENV1D && _fpUniverse->InteractionData()->RedshiftFlag()) {
      this->ComputeRedshift() ; 
      // computes _fRedshift from _fPosition : in 1D, _fPosition remains the initial position of the nu.
      _fEnergy /= (1.+_fRedshift) ;
    }
    this->Detect(aBasic, lCheck) ;
  }  
#ifdef NU_TEST
  exit(-1);
#endif
#ifdef DEBUG_OUTPUT
  std::cout<<"End of neutrino propagation" << std::endl;
#endif

  return &fParticleQueue ; 
}

vector<double> TNeutrino::CheckDetection(TBasicParam* aBasic) {
  

  vector<double> lReturn;
  lReturn.clear();
  if ( _fpUniverse->Type() == UNIVERSE_ENV1D ) {
    lReturn.push_back(1) ;
    
  } else if ( _fpUniverse->Observers()->Type() == OBSERVER_SMALLSPHERE ) {
    

    TObservers* lObs = _fpUniverse->Observers() ;
    double lDistObs;
    vector<TVector3D> lTranslators;
    TVector3D lDum ;
    TVector3D lPos=_fPosition; 
    double dt, lt; 
    TVector3D lBox; 
    TVector3D beta=_fMomentum / _fEnergy * c_light;
    TVector3D lIndex;
    lt=0;

    //Check on which side the neutrino hits the box.
    lBox.setX( (_fMomentum / _fEnergy).x() < 0? _fpUniverse->Xmin() : _fpUniverse->Xmax() );
    lBox.setY( (_fMomentum / _fEnergy).y() < 0? _fpUniverse->Ymin() : _fpUniverse->Ymax() );
    lBox.setZ( (_fMomentum / _fEnergy).z() < 0? _fpUniverse->Zmin() : _fpUniverse->Zmax() );

    lDum.set(0,0,0);
    lTranslators.push_back( lDum);
    lIndex.set(0, 0, 0); 
    int counter=0;


    while(lt < aBasic->Tmax() - _fInitTime){
#ifdef NU_TEST
      if(++counter > 10) break;
      std::cout << "lt " << lt / Mpc * c_light << "  aBasic->Tmax() " << aBasic->Tmax() / Mpc * c_light << " _fInitTime " << _fInitTime / Mpc * c_light << std::endl;
      std::cout << " aBasic->Tmax() - _fInitTime " <<  (aBasic->Tmax() - _fInitTime)/Mpc * c_light << std::endl;
#endif
      dt=           (lBox.x()-lPos.x())/beta.x() / c_light;
      dt= min(dt, ( (lBox.y()-lPos.y())/beta.y() /c_light));
      dt= min(dt, ( (lBox.z()-lPos.z())/beta.z() /c_light));
#ifdef NU_TEST
      std::cout << "beta " << beta << " dt " << dt / Mpc * c_light << "  (lBox.z()-lPos.z())/beta.z() " <<  (lBox.z()-lPos.z())/beta.z() / Mpc * c_light  << std::endl;
      std::cout << " lBox.z() " << lBox.z()/Mpc << " lPos.z() " << lPos.z()/Mpc << std::endl;
      std::cout << " lPos " << lPos / Mpc << std::endl;
#endif

      if(dt < 1e-4 *Mpc / c_light) dt=1e-4 * Mpc / c_light;
	 lPos+=dt * beta * c_light * (1 + 1e-5);

#ifdef NU_TEST
      std::cout << " lPos " << lPos / Mpc << std::endl;
#endif
      lt+=dt;

      int lflag=0; 
      if(lPos.x() > _fpUniverse->Xmax()) {
	lPos.setX(lPos.x() -  (_fpUniverse->Xmax()-_fpUniverse->Xmin())) ;
	lIndex.setX(lIndex.x()+1);
	lflag++	;
      } else if(lPos.x() < _fpUniverse->Xmin()) {
	lPos.setX(lPos.x() +  (_fpUniverse->Xmax()-_fpUniverse->Xmin())) ;
	lIndex.setX(lIndex.x()- 1);
	lflag++	;
      }

      if(lPos.y() > _fpUniverse->Ymax()) {
	lPos.setY(lPos.y() -  (_fpUniverse->Ymax()-_fpUniverse->Ymin())) ;
	lIndex.setY(lIndex.y()+1);
	lflag++	;
      } else if(lPos.y() < _fpUniverse->Ymin()) {
	lPos.setY(lPos.y() +  (_fpUniverse->Ymax()-_fpUniverse->Ymin())) ;
	lIndex.setY(lIndex.y()- 1);
	lflag++;	
      }

      if(lPos.z() > _fpUniverse->Zmax()) {
	lPos.setZ(lPos.z() -  (_fpUniverse->Zmax()-_fpUniverse->Zmin())) ;
	lIndex.setZ(lIndex.z()+1);
	lflag++;	
      } else if(lPos.z() < _fpUniverse->Zmin()) {
	lPos.setZ(lPos.z() +  (_fpUniverse->Zmax()-_fpUniverse->Zmin())) ;
	lIndex.setZ(lIndex.z()- 1);
	lflag++; 	
      }


      lDum.setX(lIndex.x() * (_fpUniverse->Xmax()-_fpUniverse->Xmin()));
      lDum.setY(lIndex.y() *(_fpUniverse->Ymax()-_fpUniverse->Ymin()));
      lDum.setZ(lIndex.z() *(_fpUniverse->Zmax()-_fpUniverse->Zmin()));

      long storA, storB;
      try{
	storA=lTranslators.size();
	storB=lTranslators.max_size();
	if(lt < aBasic->Tmax() - _fInitTime) lTranslators.push_back(lDum) ;
      }    catch (exception& e)
	{
	  std::cout << " 4rd blck " <<  e.what() << std::endl;
	  std::cout << " lTranslators.size() " << lTranslators.size() << " lTranslators.max_size() " << lTranslators.max_size() << std::endl;
	  std::cout << " lTranslators.max_size() - lTranslators.size() " << lTranslators.max_size() - lTranslators.size() << std::endl; 
	  std::cout << " - storA=lTranslators.size(); + storB=lTranslators.max_size(); " << storB - storA << std::endl;
	  exit(-1);
	}  

    }



    //    if (NEUTRINO_DETECTOR_CROWNS) {
  //     for (int ax=-1*NEUTRINO_DETECTOR_CROWNS; ax<=NEUTRINO_DETECTOR_CROWNS; ax++) {
// 	for (int ay=-1*NEUTRINO_DETECTOR_CROWNS; ay<=NEUTRINO_DETECTOR_CROWNS; ay++) {
// 	  for (int az=-1*NEUTRINO_DETECTOR_CROWNS; az<=NEUTRINO_DETECTOR_CROWNS; az++) {
// 	    lDum.set(ax*(_fpUniverse->Xmax()-_fpUniverse->Xmin()),
// 		     ay*(_fpUniverse->Ymax()-_fpUniverse->Ymin()),
// 		     az*(_fpUniverse->Zmax()-_fpUniverse->Zmin())) ;
// 	    lTranslators.push_back(lDum) ;
// 	  }
// 	}
//       }
    //   } else {
    //      lDum.set(0,0,0) ;
    //      lTranslators.push_back(lDum) ;
    //  }
#ifdef NU_TEST
    std::cout << "lt " << lt /Mpc * c_light << std::endl;
    std::cout << "lTranslators.size() : " << lTranslators.size() << std::endl;
#endif

    for (unsigned int k_trans=0; k_trans<lTranslators.size(); k_trans++) {
#ifdef NU_TEST
      std::cout << lTranslators[k_trans].x()/ Mpc << " " << lTranslators[k_trans].y()  / Mpc << " " << lTranslators[k_trans].z() / Mpc << std::endl;
#endif

      for (int j=0; j<lObs->Nb(); j++) {
	//	std::cout << lObs->Positions(j) / Mpc << std::endl;
	TVector3D lObsPart = _fInitPosition - (lObs->Positions(j)+lTranslators[k_trans]) ;
	lDistObs = lObsPart.mag() - lObs->Radius() ;
	if (lDistObs<=0) lReturn.push_back(0);
	else {
	  double lXP, lDiscrim ;
	  lXP = lObsPart.dot(_fMomentum.unit()) ;
	  lDiscrim = pow(lXP,2) - lObsPart.mag2() + pow(lObs->Radius(),2) ;
	  if ( lDiscrim >= 0 && lXP <= 0) { 
	    // Direction defined by momentum hits sphere: neutrino hits the detector.
	    lReturn.push_back(-lXP); //= 1 ;
	    //_fPosition -= _fMomentum.unit()*lXP ; // Important correction to position
	  }
	}
      }
    } 
  } else if ( _fpUniverse->Observers()->Type() == OBSERVER_LARGESPHERE ) {
    TObservers* lObs = _fpUniverse->Observers() ;
    double lDistObs ;
    int flag=0;
    int j=0 ; // Only one observer allowed
    lDistObs = (_fPosition-_fSourcePosition).mag() - lObs->Radii(j) ;
    if (lDistObs >= 0) flag = 0; // particle out of the sphere
    else flag = 1 ; // particle inside the sphere
    if (flag) {
      double lCosTheta = ((_fPosition-_fSourcePosition).unit()).dot(_fMomentum.unit()) ;
      double lSP = (_fPosition-_fSourcePosition).mag() ;
      double lDist = sqrt(pow(lObs->Radii(j),2) - lSP*lSP*(1-lCosTheta*lCosTheta)) - lSP*lCosTheta ;
      lReturn.push_back(lDist);
      //      _fPosition += _fMomentum.unit() * lDist ;
    }
    
  } else if ( _fpUniverse->Observers()->Type() == OBSERVER_NO ) {
    // lReturn = 0 ;
    lReturn.clear();
    
  } else throw TCrpErr("Error : CheckDetection for photons not implemented in this configuration" );
  

  return lReturn ;


}

void TNeutrino::Detect(TBasicParam* aBasic, vector<double>lCheck) const {
  for(int i=0; i < lCheck.size(); i++){
  TOutputData* lOutput = aBasic->OutputData() ;
  if (_fpUniverse->Type() == UNIVERSE_ENV1D)
    lOutput->Add1DNeutrino(_fFlavor,_fSourcePosition.mag(),_fInitPosition.mag(),_fInitialEnergy,_fEnergy,_fSourceEnergy, _fInitType) ;
  else if (_fpUniverse->Type() == UNIVERSE_LARGESCALE){
    if(lCheck[i] < (aBasic->Tmax() - _fInitTime)*c_light)
      lOutput->Add3DNeutrino(_fFlavor,_fSourcePosition,_fInitPosition,_fPosition + _fMomentum.unit() * lCheck[i],_fMomentum,_fSourceEnergy, _fInitType) ;}  
  else throw TCrpErr( "Detection not implemented in this environment" );
  }
}

