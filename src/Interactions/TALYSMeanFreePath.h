#ifndef _TALYSMeanFreePath_H_
#define _TALYSMeanFreePath_H_

#include <cmath>
#include <vector>
#include <TabulatedTALYSY.h>
#include <iostream>
#include <iomanip>
#include <CMB.h>
#include <IR.h>
#include "nucleusdb.h"

#include <TGraph.h>
#include <TCanvas.h>

#include "CLHEP/Vector/ThreeVector.h"
#include <vector3d.h>

/**
   @class TALYSMeanFreePath

   @author  Nils Nierstenhoefer nierstenhoefer@physik.uni-wuppertal.de

   @brief Prototype class to calculate and provide the mean free path for photo disintegration of UHE nuclei. 

   <strong> The General Organization of the Photo Disintegration </strong>
   
   The probability P for a reaction to take place after a distance 
   \f$\Delta \f$x can be sampled from an exponential distribution.
   using the mean free path \f$\lambda\f$ 

   \f$P=1-\exp (-\frac{\Delta x }{ \lambda}). \f$ 

   \f$\lambda\f$ for a photo disintegration reaction of a nucleus with an ambient photon density 
   \f$n(\epsilon, z)\f$ can be found by calculating 
   
    \f$
    \lambda^{-1} 

     =

    \frac{1}{2 \gamma^2 }
    \int_{0}^{ \infty } \frac{n(\epsilon , z)}{\epsilon^2}
    
    \left\{ 
    \int_{0}^{2 \gamma \epsilon} \epsilon ' \sigma (\epsilon ' )\,d\epsilon ' 
    \right\} 

    d\epsilon
    \\
     =

    \frac{1}{2 \gamma^2 }
    \int_{0}^{ \infty } \frac{n(\epsilon , z)}{\epsilon^2}
    \cdot
    \sigma_{\rm{avrg}}( \epsilon')
    \cdot
    d\epsilon.

    \label{AttLengthInt}
    \f$

    Here, \f$\epsilon\f$ is the photon energy in the laboratory frame and \f$\gamma\f$ is the Lorentz factor 
    of the nucleus. The prime marks values in the nucleus' rest frame. \f$\sigma \f$ 
    is the cross section (it is assumed that the TALYS cross sections are in given in the nucleus' rest frame).

    The photo disintegration becomes dominant if the \f$\epsilon '\f$ is of the order
    the binding energy ~8-9 MeV/nucleon. Although, the main photo disintegration channels are the loss of one proton or 
    neutron, at higher energies channels with larger multiplicities might be excited. This includes the loss of 
    heavier particles. Here, the TALYS framework was used to calculate the photo disintegration cross sections
    of nuclei with mass numbers A>6 and neutron number N>2 in the energy range of \f$[10^{-3}, 250]~\rm{MeV}\f$. This includes the 
    combined losses of
    proton (p), neutron (n), deuterium (d), tritium (t), helium-3 (He-3) and \f$\alpha\f$-particles. 
    For A<10 the TALYS predictions become unreliable and the cross sections of the most stable nuclei 
    were replaced with the ones according to Rachen et al., Geant-4 or available experimental data.     

    Independently of the application, the inner integral \f$\sigma_{\rm{avrg}}\f$ can be calculated before a 
    simulations is started, and it can be tabulated as a function of the photon energy \f$\epsilon'=2\gamma \epsilon\f$. As it merely 
    is an averaging of the cross section for an mono-energetic, 
    isotropic photon density in the laboratory frame, \f$\sigma_{\rm{avrg}}\f$ will be referred to as 
    "averaged cross section" later on - although the units are not m^2. In case of variable photon fields the outer integral
    can then be calculated on the fly.  
    In case of a constant photon density even the mean free path \f$\lambda\f$ itself can be 
    tabulated as function of the Lorentz factor \f$\gamma\f$ of the nucleus.    
   
    There are two prototype classes which are linked with the calculation and handling
    of the mean free path \f$\lambda\f$:   
    <ul>
    <li> TabulatedTALYSY
    is the prototype class for the handling of the photo disintegration cross section \f$\sigma (\epsilon' )\f$, 
    the average cross section \f$\sigma_{\rm{avrg}}(\epsilon'=2\gamma \epsilon )\f$ or the inverse mean 
    free path \f$\lambda(\gamma )^{-1}\f$ as "tabulated" functions using the derived classes
        <ol>
        <li> TabulatedTALYSY:TabulatedTALYSCrossSection,
        <li> TabulatedTALYSY:TabulatedTALYSAveragedCrossSection,
        <li> TabulatedTALYSY:TabulatedTALYSMeanFreePath,
	  <ul>
	     <li> TabulatedTALYSMeanFreePath:TabulatedTALYSMeanFreePath_CMB,
	     <li> TabulatedTALYSMeanFreePath:TabulatedTALYSMeanFreePath_IRB,
	     </ul>
        </ol>
	respectively.
     <li> TALYSMeanFreePath 
     is the prototype class to calculate or provide \f$\lambda(\gamma )\f$ using the data classes 
     derived from TabulatedTALYSY. For this purpose, it calculates \f$\lambda(\gamma )\f$
     by solving eq. (#2) with \f$\sigma (\epsilon' )\f$ from TabulatedTALYSCrossSection,
     by solving only the outer integral  eq. (#3) using \f$\sigma_{\rm{avrg}}(\epsilon')\f$ from
     TabulatedTALYSAveragedCrossSection or immediately returns \f$\lambda(\gamma )\f$ from 
     TabulatedTALYSMeanFreePath using the derived classes
        <ol>     
        <li> TALYSMeanFreePath:TALYSMeanFreePathAccurate,
        <li> TALYSMeanFreePath:TALYSMeanFreePathAvrgd,
        <li> TALYSMeanFreePath:TALYSMeanFreePathTabulated,
        </ol>
	respectively.
     </ul>


The Monte Carlo procedure itself is implemented in TNucleus and TSophiaInteractions via the functions
<ol>     
    <li> TSophiaInteractions::GetPDTimeStep (...),
    <li> TSophiaInteractions::GetExclusivePDChannelDirectly(...),
    <li> TNucleus::TALYSNucleiPhotoDisintegration(...).
</ol>   

The xml-configuration for the photodisinetgration is done in TSophiaInteractions::TSophiaInteractions(). 

<strong> Labeling of Nuclei and Reaction Channels </strong>

Details are explained in TabulatedTALYSY.

<strong> Comments on the TALYSMeanFreePath Prototype Class </strong>

As already stated above, this is a prototype class to provide the mean free path
for photo disintegration reactions of nuclei. 

It turned out to be useful to first perform
a substitution \f$t= \ln(\epsilon/\epsilon_{0}) \f$ in the outer integral of equation (2.) before
calculating it. 
Furthermore, the integral is solved with an 8th order Gauß-Legendre algorithm (8 sample points). 
But, higher precision is needed to guarantee “overall” convergence for all nuclei. Thus,  
the integration interval is disjoint into \f$N\f$ pieces which become integrated separately(linearity of integral).   
   
\f$
     \lambda^{-1} 
     
     = 
     
     \frac{1}{2 \gamma^2 }
     \int_{t_{\rm{min}}}^{t_{\rm{max}}} 
        \frac{ n(\epsilon (t) , z)}{\exp^2(t)\epsilon_{0}^2}
        \cdot
        \sigma_{\rm{avrg}}( \epsilon' (t))
        \cdot
	\exp(t)\epsilon_{0}
	dt

     \\

     = 

     \frac{1}{2 \gamma^2}
     \sum_{i=1}^{N}
     \int^{t_{\rm{min}}+(t_{\rm{max}}-t_{\rm{min}})i/N}
     _{t_{\rm{min}}+(t_{\rm{max}}-t_{\rm{min}})(i-1)/N} 
     \frac{n(\epsilon (t) , z)}{\exp(t) \epsilon_{0}}
        \cdot
        \sigma_{\rm{avrg}}( \epsilon' (t))
        \cdot
	dt
\f$
     
Here, \f$\epsilon(t)=\exp(t)\epsilon_{0}\f$ and 
\f$t_{\rm{min}/\rm{max}}=\ln(\epsilon_{\rm{min/\rm{max}}}/\epsilon_{0})\f$ 
where \f$\epsilon_{\min}\f$ and \f$\epsilon_{\max}\f$ 
are separately fixed for the CMB \f$(\epsilon_{\min}=4\,10^{-19}\f$ GeV, \f$\epsilon_{\max}=10^{-11}\f$ GeV) and 
IRB \f$(\epsilon_{\min}=10^{-12}\f$ GeV, \f$\epsilon_{\max}=10^{-7}\f$ GeV) for the mean free path table 
calculations. In principle the photo disintegration routines can handle 
variable photon fields by using the tabulated averaged cross sections and solving the outer integral 
on the last equation on the fly. But currently this is only an undocummented feature in this version of CRPropa.
In spite of that the redshift dependency is approximately taken into account by a scaling function \f$s(z)=(1+z)^3\,\,e(z)\f$.
The functional for of \f\f$e(z)=1\f$ in case of the IRB is discussed in the CRPropa 2.0 (or in the doxygene: IRBzEvolutionFactory) 
paper while \f$e(z)=1\f$ in the case of the CMB. Using \f$s(z)\f$, the mean free path scales as

\f$ \lambda^{-1}=(1+z)^3 \cdot \lambda[(1+z)E,z=0]^{-1}\f$.

combing the last two equations and using \f$\gamma=E/M\f$ one finds 

\f$
     \lambda'^{-1} 
     =
     \frac{M^2(1+z)^3}{2 E^2(1+z)^2}\,e(z)
     \sum_{i=1}^{N}
     \int^{t_{\rm{min}}+(t_{\rm{max}}-t_{\rm{min}})i/N}
     _{t_{\rm{min}}+(t_{\rm{max}}-t_{\rm{min}})(i-1)/N} 
     \frac{n[\exp(t) \epsilon_{0}(1+z) , z]}{\exp(t) \epsilon_{0}}
        \cdot
        \sigma_{\rm{avrg}}[ 2\cdot  \exp(t) \epsilon_{0} (1+z)\epsilon /M ]
        \cdot
	dt
	\label{eq:MFPFinal}
\f$

This is exactly the form which is used in the code - that justifies the lenghty notation.

<strong> Units </strong>

Assuming \f$\hbar=c=1\f$, the units of equation 2 can be found to be

\f$
[\lambda']^{-1}=  1 \cdot \frac{\rm{GeV}^2}{\rm{GeV}^2} \cdot \left\{ \rm{GeV} \cdot mbarn \cdot \rm{GeV}\right\} 
            \cdot\rm{GeV} 
         =  \rm{GeV}^3 \cdot mbarn
\f$

where the units of the photon fields [n]=\f$\rm{GeV}^2\f$ and energies are measured in units of GeV. \f$\hbar=1\f$
and \f$c=1\f$ can be written as m=s/c' and s=1/6.582e-25/GeV where c' is the speed of light (without units). Hence, 
\f$ \rm{mbarn}=10^{-31}\rm{m}^2=10^{-31}(\rm{s}/\rm{c'})^2=2.568\cdot \rm{GeV}^{-2} \f$.

\f$
[\lambda']^{-1}=2.568\cdot \rm{GeV}
\f$
 
This can be rewritten as SI length in Mpc: \f$\rm{GeV}=\sqrt{2.568/10^{-31}}\rm{m}^{-1}
=5.068\cdot 10^{15} \cdot 3.068\cdot 10^{22} \rm{Mpc}^{-1}\f$ 

\f$
[\lambda']^{-1} =  2.568\cdot 5.068\cdot 10^{15} \cdot 3.068\cdot 10^{22} \rm{Mpc}^{-1} \\
                =  1.301\cdot 10^{16}  \cdot 3.068\cdot 10^{22} \rm{Mpc}^{-1}
\f$

In other words, one can solve the integral of equation 2 using \f$\hbar=c=1\f$, energies and masses
in GeV, the cross section in mbarn and the photon density in units of \f$\rm{GeV}^2\f$. In a 
second step it can be converted to the wanted unit - \f$\rm{Mpc}^{-1}\f$ - by multiplying
it with the factor \f$D=1.301\cdot 10^{16}  \cdot 3.068\cdot 10^{22} \rm{Mpc}^{-1} \rm{GeV}^{-1}\f$ 

<strong> Plotting the Mean Free Path or the Cross Section for Photo Disintegration </strong>

Up to now there is no xml-tag which will cause CRPropa to create plots of the mean free path or cross section for 
photo disintegration as function 
of the nucleus' energy or gamma factor. But, the needed lines of code can "be commented in by hand" in the construtor 
TNucleus::TNucleus() in nucleus.cc by adding
CPPFLAGS="-DPLOT_MEAN_FREE_PATH_ROOT"
to the configure call. After recompilation, CRPropa will then create root files with the mean
free path or cross section plots if it is started in a xml configuration which is linked with propagation of nuclei - 
this will call the constructor. A set of root files will be created.

<strong> Todos and Notes </strong>
<ul>
<li> Functions that are handed to the gauss( . . . ) integration routine have always the same set of arguments and sometimes not all of them are needed. 
More precisely: functs_int_times_PhtnBckGrnd(. . .) has unneeded arguments such that it can be given to gauss(. . . ) as an argument.
<li> Currently the photo disintegration (PD) routines create their own instances of the photon fields in variable IRB mode (in sophiainteractions.cc).
But the IRB is already constructed while creating the TUniverse and ,indeed , this instance of the IRB should be handed to the PD code. E.g. give the photonfields via argument to TSophiaInteractions::GetExclusivePDChannelDirectly(...) to the PD routines via TNucleus::TALYSNucleiPhotoDisintegration(...) - note, a pointer to the TUniverse exits in TNucleus, but not in TSophiainteractions.
<li> All the PD classes could inherit from TXmlParam. Then, one could remove the xml configuration from the 
constructur of sophiainteractions and move it to the corresponding PD classes' constructors.   
<li> Get rid of the CPPFLAGS "-DPLOT_MEAN_FREE_PATH_ROOT" flag to compile CRPropa to include the plotting routines for PD - introduce a corresponing XML keyword instead.
<li> Implementation of a higher order numerical integration.
<li> \f$\epsilon_{0}\f$, \f$t_{\rm{min}}\f$  and \f$t_{\rm{max}}\f$ should be given as a 
parameter to the member call of Total_rate_log(...). Note, in the moment \f$\epsilon_{0}\f$
is provided by the CMB and IR classes while \f$t_{\rm{min}}\f$ and \f$t_{\rm{max}}\f$ are
hard coded and have the same value for the case of CMB and IR.
</ul>
*/


class TALYSMeanFreePath {

public:
  
  typedef double (TALYSMeanFreePath::*ClassMemPtr) (double,
						    TVector3D,
						    unsigned long int,
						    unsigned long int,
						    double,
						    double,
						    double,
						    int);
TALYSMeanFreePath();
//This constructor is for the case of the averaged cross section tables. For the on the fly folding with the photonfields photonenergy range (epsilon_min, epsilon_max) and epsilon_zero are needed. Note, epsilon is for the substitution of t=ln(e/epsilon_zero) where epsilon is the photon energy in the lab frame - c.f. doxygene page.
TALYSMeanFreePath(double epsilon_zero, double epsilon_min, double epsilon_max);
~TALYSMeanFreePath();

  virtual void SetNIntegratorCalls(int NIntegrationCalls);
  /**< Number of joints (= number of gauss integrator calls) to solve outer integral to calculate mean free path (cf N in equation 1 of section "Comments on the TALYSMeanFreePath Prototype Class" above). */

  virtual double functs(double s, 
			TVector3D Position,
			unsigned long int InputNucleusId,
			unsigned long int OutputId,
			double NucleusMass,
			double EnergyOfNucleusGeV,
			double z,
			int TotalOrExclusive)=0;
  /**< Returns \f$\epsilon'\sigma(\epsilon')\f$. If TotalOrExclusive=0 calculate total channel. If TotalOrExclusive=1 calculate exclusive channel.*/

  
  virtual double functs_int(double UpperBorder,
			    TVector3D Position,
			    unsigned long int InputNucleusId,
			    unsigned long int OutPutNucleusId,
			    double NucleusMass,
			    double EnergyOfNucleusGeV,
			    double z,
			    int TotalOrExclusive);
  /**< Returns \f$\int \epsilon'\cdot \sigma(\epsilon')\f$. If TotalOrExclusive=0 calculate total channel. If TotalOrExclusive=1 calculate exclusive channel. Some arguments are nor needed. They were added to meet the arguments list for the gauss() routine.*/


  virtual double Total_rate(std::vector<PhotonBackground*> PhotonBkgrd,
			    std::vector<TabulatedTALYSY*> tTabulatedTALYSY_Vec,
			    PD_TabDataType tPD_TabDataType,
			    TVector3D Position,
			    unsigned long int InputNucleusId,
			    unsigned long int OutPutNucleusId,
			    double NucleusMass,
			    double EnergyOfNucleusGeV,
			    double z,
			    int TotalOrExclusive);
    /**< Returns \f$\lambda^{-1}\f$ with substitution \f$\epsilon = \epsilon_{0}(1-t)/t\f$ in outer integral. A vector of photon fields can be given. It should be a vector with no elements in case of the mean free path tables. If TotalOrExclusive=0 calculate total channel. If TotalOrExclusive=1 calculate exclusive channel.*/

  virtual double functs_int_times_PhtnBckGrnd(double t,
					      TVector3D Position,
					      unsigned long int InputNucleusId,
					      unsigned long int OutPutNucleusId,
					      double NucleusMass,
					      double EnergyOfNucleusGeV,
					      double z,
					      int TotalOrExclusive);
  /**< Returns \f$ \int n(z, \epsilon) / \epsilon^2 * \int_{0}^{2 \gamma \epsilon} \epsilon'\cdot \sigma(\epsilon')\f$ where \f$n(z, \epsilon)\f$ is the photon density. This is the case of the substitution \f$\epsilon = \epsilon_{0}(1-t)/t\f$. If TotalOrExclusive=0 calculate total channel. If TotalOrExclusive=1 calculate exclusive channel.Some arguments are nor needed. They were added to meet the arguments list for the gauss() routine.*/
  
   
  virtual double Total_rate_log(std::vector<PhotonBackground*> PhotonBkgrd,
				std::vector<TabulatedTALYSY*> tTabulatedTALYSY_Vec,
				PD_TabDataType tPD_TabDataType,
				TVector3D Position,
				unsigned long int InputNucleusId,
				unsigned long int OutPutNucleusId,
				double NucleusMass,
				double EnergyOfNucleusGeV,
				double z,
				int TotalOrExclusive);
  /**< Returns \f$\lambda^{-1}\f$ with substitution \f$t=ln(\epsilon / \epsilon_{0})\f$ in outer integral. A vector of photon fields can be given. It should be a vector with no elements in case of the mean free path tables. If TotalOrExclusive=0 calculate total channel. If TotalOrExclusive=1 calculate exclusive channel.*/

  virtual double functs_int_times_PhtnBckGrnd_log(double t,
						  TVector3D Position,
						  unsigned long int InputNucleusId,
						  unsigned long int OutPutNucleusId,
						  double NucleusMass,
						  double EnergyOfNucleusGeV,
						  double z,
						  int TotalOrExclusive);
  /**< Returns \f$ \int n(z, \epsilon)/ \epsilon^2 * \int_{0}^{2 \gamma \epsilon} \epsilon'\cdot \sigma(\epsilon')\f$ where \f$n(z, \epsilon)\f$ is the photon density. This is the case of the substitution \f$t=ln(\epsilon / \epsilon_{0})\f$. If TotalOrExclusive=0 calculate total channel. If TotalOrExclusive=1 calculate exclusive channel. Some arguments are nor needed. They were added to meet the arguments list for the gauss() routine.*/
  


  virtual double MFPNanInfOrSmallerZeroCheck(double MFP);
  /**< Returns std::numeric_limits<double>::min() if MFP is nan, inf or zero. */
  
  void RootFileWithAllMFPs(std::vector<PhotonBackground*> PhotonBkgrd,
			   std::vector<TabulatedTALYSY*> tTabulatedTALYSY_Vec,
			   TVector3D Position,
			   double fXMin,
			   double fXMax,
			   double fXNSteps,
			   double redshift);
  /**< Creates a set of files which hold plots of \f$\lambda^{-1}\f$ for all channels.*/
  
  TGraph* RootPlotMFPVsX(std::vector<PhotonBackground*> PhotonBkgrd,
			 std::vector<TabulatedTALYSY*> tTabulatedTALYSY_Vec,
			 TVector3D Position,
			 unsigned long int InputNucleusId,
			 long unsigned int OutId,
			 double fXMin,
			 double fXMax,
			 double fXNSteps,
			 double redshift,
			 int TotalSumOrExclusive);
  /**< Plotting routine.*/
  
  TGraph* RootPlotExclSumMFPVsX(std::vector<PhotonBackground*> PhotonBkgrd,
				std::vector<TabulatedTALYSY*> tTabulatedTALYSY_Vec,
				TVector3D Position,
				unsigned long int InputNucleusId,
				double fXMin,
				double fXMax,
				double fXNSteps,
				double redshift);
  

  int  _fNIntegrationCalls;
  double _fEpsilon0;
  double _fEpsilon_min;
  double _fEpsilon_max;
  PhotonBackground* fPhotonBkgrd;

 protected:
  virtual std::vector<unsigned long int> GetCombinedExclusiveChannelList(std::vector<TabulatedTALYSY*> tTabulatedTALYSY_Vec,
						 unsigned long int NucleusId);
  virtual std::vector<unsigned long int> GetUniqueVecFrom2Vec(std::vector<unsigned long int> vec1, std::vector<unsigned long int> vec2);

 template < typename T > double gauss (T function,
					TVector3D Position,
					double A,    
					double B,
					unsigned long int InputNucleusId,
					unsigned long int OutPutNucleusId,
					double NucleusMass,
					double EnergyOfNucleusGeV,
					double z,
					int TotalOrExclusive){
  
  if(B<A){
    std::cout<<"Warning B<=A in GAUSS."<<std::endl;
    return(0);
  }

  std::vector<double> Nodes;
  Nodes.push_back(0.0950125098);
  Nodes.push_back(0.2816035507);
  Nodes.push_back(0.4580167776);
  Nodes.push_back(0.6178762444);
  Nodes.push_back(0.7554044083);
  Nodes.push_back(0.8656312023);
  Nodes.push_back(0.9445750230);
  Nodes.push_back(0.9894009349);

  std::vector<double> Weights;
  Weights.push_back(0.1894506104);
  Weights.push_back(0.1826034150);
  Weights.push_back(0.1691565193);
  Weights.push_back(0.1495959888);
  Weights.push_back(0.1246289712);
  Weights.push_back(0.0951585116);
  Weights.push_back(0.0622535239); 
  Weights.push_back(0.0271524594);

  double XM = 0.5*(B+A);
  double XR = 0.5*(B-A);
  double SS = 0.;

  //TODO: Take the method from, "Numerical methods using MathLab", 4th edition, Mathews, Fink which needs just one function call -> faster!?
  for(int NJ=1;NJ<=8;NJ++){
    double DX = XR*Nodes[NJ-1];

    SS = SS + Weights[NJ-1] * ((this->*function)(XM+DX, Position ,InputNucleusId, OutPutNucleusId, NucleusMass, EnergyOfNucleusGeV, z, TotalOrExclusive) 
			       + 
			       (this->*function)(XM-DX, Position ,InputNucleusId, OutPutNucleusId, NucleusMass, EnergyOfNucleusGeV, z, TotalOrExclusive));	
  }
  
  return(XR*SS);  
 }
  
}; 
#endif
