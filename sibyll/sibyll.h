/** 
	@file     sibyll.h
        @brief    Definition file for fortran / C interface
*/

#ifndef _SIBYLL_H_
#define _SIBYLL_H_

extern "C" {
  void sibyllevent_(int&, double&, double[][15000], int[], int&) ;
}

/*
The arguments are the following (see sibyll_interface.f) :

     INPUT
     iproj = 0 -> p ; 1 -> n
     Ein : in GeV (SIBYLL standard energy unit)
   
     OUTPUT
     OutPart,OutPartType,NbOutPart = output data:
     P(8000,5) list of 4-momenta + masses of output particles
     LList(8000) list of output particle IDs
     NP nb of output particles
    

cc      13  proton
cc      14  neutron
cc      -13 antiproton
cc      -14 antineutron
cc      1   photon
cc      2   e+
cc      3   e-
cc      15  nu_e
cc      16  antinu_e
cc      17  nu_muon
cc      18  antinu_muon

Warning with passing arguments (variables->pointers) in fortran->C !

*/

#define sibyllevent sibyllevent_

#endif
