c****************************************************************************
c
c   SIBYLLEVENT
c
c   interface between Sybill and CRPropa
c  
c****************************************************************************

      subroutine sibyllevent(iproj,Ein,OutSib,OutSibType,NbOutSib)
      
c*****************************************************************************
c     INPUT
c     iproj = 0 -> p ; 1 -> n
c     Ein : in GeV (SIBYLL standard energy unit)
c     
c     OUTPUT
c     OutPart,OutPartType,NbOutPart = output data:
c     P(15000,5) list of 4-momenta + masses of output particles
c     LList(15000) list of output particle IDs
c     NP nb of output particles
c    
c*****************************************************************************


    
      implicit double precision(A-H,O-Z)
      
c      integer Np,List(2000)
c      real P(2000,5)


      common/S_PLIST_SIB/P(15000,5),List(15000),Np 
      common/S_MASS1/ AM(49), AM2(49)
      common/decay/i_decay
      
      double precision Ein,gamma, beta, pm, pm_targ
      
      double precision OutSib(15000,5)
      integer OutSibType(15000),NbOutSib

      if (iproj.eq.0) then 
         L0=13 
      else if (iproj.eq.1) then
         L0=14
      else
         write(*,*) 'sibyllevent: incorrect incoming particle'
         stop
      endif
      
      call sibyll_ini           ! Initialization of cross sections

c     Sibyll parameters
c     ------------------
      i_decay = 1               ! 1:decay of secondaries ON, 0:decay OFF
      
c     mass number A of the target nucleus 
      iatarg = 1                ! proton
c     iatarg=0 is an "air" nucleus  (superposition of oxygen and nitrogen)
      
      pm = AM(L0)
      pm_targ = AM(13)          ! proton
      
c     Kinematic variables
c     -------------------
      Ecm = sqrt(2d0*pm_targ*Ein + pm**2 + pm_targ**2) ! Ecm = sqrt(s)
      gamma = Ecm/(2d0*pm) 
      beta = sqrt(1d0-1d0/gamma**2)
      
      call sibyll(L0,iatarg,Ecm)
      
      icont = 0
      do i=1,Np
         if(abs(List(i)).LT.10000) then
            icont = icont + 1
            do j=1,5
               OutSib(icont,j)=P(i,j)
            enddo
c     Lorentz transformation
            OutSib(icont,4) = gamma*P(i,4) + gamma*beta*P(i,3)
            OutSib(icont,3) = gamma*P(i,3) + gamma*beta*P(i,4)
            OutSibType(icont)=List(i)
         else 
            cycle
         endif
      enddo
      NbOutSib=icont

      return
      end
