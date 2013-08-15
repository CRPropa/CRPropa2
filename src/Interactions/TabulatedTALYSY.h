#ifndef _TabulatedTALYSY_H_
#define _TabulatedTALYSY_H_

#include <fstream>
#include<iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <nucleusdb.h>
#include <sstream>

#include <TGraph.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TLegend.h>
#include <TH1F.h>
#include <TFile.h>
#include <TGraph.h>
#include <PhotonBackground.h>
#include <vector3d.h>

/**
   @class TabulatedTALYSY

   @author  Nils Nierstenhoefer nierstenhoefer@physik.uni-wuppertal.de

   @brief Prototype class to handle data needed for photo disintegration of UHE nuclei. 

   <strong> The General Organization of the Photo Disintegration </strong>
   
   The details are explained in the documentation of TALYSMeanFreePath.

   <strong>Labeling of Nuclei and Reaction Channels </strong>

Nuclei are labeled with an id \f$I(A,Z)\f$ according to the following convention within CRPropa
\f$I(A,Z)=Z*1000+A\f$
where Z, A are the charge, mass number.

A nuclear transition  \f$I(A,Z)\rightarrow I(A',Z')+ \sum_{i} I(A_{i},Z_{i})\f$ is linked 
with a change in the charge- and mass number. Here, the \f$I_{i}\f$ are the set of nuclei
which were disintegrated from the initial nucleus. Let's define \f$\Delta A = A'-A\f$ and
\f$\Delta Z = Z'-Z\f$. All possible channels for the aforementioned transition are combination 
\f$(\#n,\#p,\#d\#t,\#He-3,\# \alpha )\f$ which obey \f$\sum_{i} A_{i}=\Delta A\f$ 
and \f$\sum_{i} Z_{i}=\Delta Z.\f$. \#n,\#p,\#d, etc. means number of neutrons, protons, deuterium 
and so forth. These reaction channels are called "exclusive" channels. E.g. the two 
exclusive channels (2,2,0,0,0,0) and (0,0,0,0,0,1) would cause the same nucleus transition 
but are linked with different reaction products and probabilities. The exclusive channels are
labeled with an additional unique id

\f$ E(\#n,\#p,\#d,\#t,\#He-3,\# \alpha )=  \#n\cdot 10^6 + \#p\cdot 10^5 + . . . +\# \alpha\cdot 10^0\f$

<strong> Installation of the Photo Disintegration Cross Section Data </strong>

The cross section data package for the photo disintegration does not come with the CRPropa framework. This is mainly due to
its size (~1 GB). The package can be downloaded and installed by the use of

<CODE> ./GetPDCrossSections.sh </CODE>

This script is located in the trunk folder and will install the recommended photo disintegration tables with a thinning 
factor of 90% (where 100% means no thinning) by default. If you want to use e.g. the cross section without thinning start the
install script anyway. Then copy the file from CrossSectionPackage/TabulatedTALYS*_thn100/ into the corresponding trunk
folder

<CODE> 
cp CrossSectionPackage/TabulatedTALYSAveragedCrossSection_XSBins500_ABins200_thn100/PD*.cax TabulatedTALYSAveragedCrossSection/ \n
cp CrossSectionPackage/TabulatedTALYSMeanFreePath_XSBins500_ABins200_MFPBins200_thn100/PD*.cmt TabulatedTALYSMeanFreePath/ \n
make install
</CODE>  

You can create tables with a thinning level of your own choice by following the instructions given in the section 
"Creation of the Data Tables" (s. below). 

<CODE>make install</CODE> will copy the the data into the install directory. 

<strong> Comments on the TabulatedTALYSY Prototype Class </strong>

A continuous function \f$Y(x)\f$ can be approximated by a "tabulated" function, that is: Store pairs 
\f$(x_{i},Y(x_{i}))\f$ and fake a continuous function by performing a linear approximation between those 
pairs. This functionality is the main task of TabulatedTALYSY prototype class for the photo disintegration
data (cross section \f$\sigma\f$ , averaged cross section \f$\sigma_{\rm{avrg}}\f$ or the mean free 
path \f$\lambda\f$ itself). Due to the high multiplicity of exclusive channels, there is a huge amount 
of data and accessing it is a bottleneck.  Therefore, a simple indexing system was introduced which comes
along with the distribution of the data into a set of five data files:
<ul>     
     <li> <strong> *InitNucleusId.c*  </strong>
     
     The starting line of this file holds:
         <ol> 
	    <li>the data type identified by 0,1 or 2 which corresponds to \f$\sigma\f$, \f$\sigma_{\rm{avrg}}\f$ 
	    or \f$\lambda\f$.
	    <li> Minimum x value, 
	    <li> maximum x value,
	    <li> number of bins (= number of pairs \f$(x_{i},Y(x_{i}))\f$),
	    <li> thinning factor (only given in case of \f$\sigma_{\rm{avrg}}\f$ 
	    and \f$\lambda\f$).
	 </ol>
      The following lines hold the ids \f$I(A,Z)\f$ of the nuclei which can be tracked in CRPropa and the rows 
      where the reaction channels which are linked with the nuclei can be found in the *Excl*CrossId.c*
      and *Total*Cross.c* files. The 2nd number correspond the total- and the 3rd number to
      the exclusive reactions channels.

     <li>  <strong> *Excl*CrossId.c* </strong> 

     The rows in this files begin with the id of an exclusive channel \f$E(\#n,\#p,\#d,\#t,\#He-3,\# \alpha )\f$.
     It is followed by a pairs of numbers. The pair gives start- and end row of the \f$Y_{i}\f$ values 
     of the exclusive channel. The remaining two lines
     hold the loss in charge \f$\Delta Z\f$ and mass  \f$\Delta A\f$ which is linked with this exclusive
     channel. 

     <li>  <strong> *Excl*Cross.c* </strong> 

     This files finally holds the \f$Y_{i}\f$ values in the rows as specified in the *Excl*CrossId.c* files (
     ordered by increasing \f$x_{i}\f$) values.

     <li>  <strong> *ExclSum*Cross.c*  </strong>

     For each isotope the sum of all exclusive channels is given in this file. This data can be used to provide/calculate the total reactions rate
     in TSophiaInteractions::GetPDTimeStep (...) directly without summing over all exclusive channels during runtime ( ==> much faster especially if the folding 
     with the photon field is calculated on the fly!).

     This file is not needed and thus not given in case of the cross section \f$\sigma\f$.

</ul> 

For example, what is the cross section for an oxygen nucleus to loose one neutron and one \f$\alpha\f$ particle in a photo disintegration 
reaction with a photon of an energy of 150 MeV? Oxygen has nucleus id of I(16,8)=8016. The id of the exclusive channel is 
\f$E(1,0,0,0,0,0,1)=100001\f$. Assume one has a table with NBins=500 bins in the energy range from [0.001, 250) MeV.
Then, the bin width is ~0.5 MeV. Hence, the cross section for 150 MeV can be calculated using a linear approxmation between
the pairs \f$(X_{i}, Y_{i})\f$ with i=300 and i=301. The two pairs can be found in the data files by first locating the
line which start with "8016" in the PDInitNucleusId.cxs (xs= cross section). The line contains five numbers:

<CENTER>
<CODE>
...

8016  6374 6646

...
</CODE>
</CENTER>

As we are looking for the exclusive channel 100001 we have to deal with the last two values. There are 
6646-6374=272 available exclusive channels for oxygen (8016). They are listed in the files PDExclCrossId.cxs in 
the rows 6646 to 6374. The corresponding lines might look like that.

<CENTER>
<CODE>
111101 39279500 5 11\n
041100 39280000 6 9\n
100101 39280500 3 8

. . . 

<strong> 
100001 39408500 2 5 
</strong>

. . . 

431000 39414500 4 9\n
012100 39415000 4 8\n
101001 39415500 3 7\n
</CODE>
</CENTER>

Thus the corresponding \f$Y_{i}\f$ values needed for the linear approximations can be found in line 39408500+300 
and line 39408500+301 in the file PDExclCross.cxs. Column three and for hold the charge- and mass loss linked
with the exclusive channel. 

Note: for \f$\sigma_{\rm{avrg}}\f$ and \f$\lambda\f$ there are two additional data columns. They hold the first and last line of 
where the total cross section for the isotope can be found in the corresponding *ExclSum*Cross.c* file.

<strong>Singelton classes</strong>

This class is organized as a "singleton". That means one can only create one instance
 of this class which is accesible from all over the programm. This is clearly useful due to the cross-section 
data which has to be red into the memory only once and only if needed. It follows a short sketch which shows the structure of a class which has this 
feature (extracted from wikipedia):   

<CODE>
class N \n
 {
 public: \n
    \t static N* instance() \n
    \t{ \n
       \t\t if(! _instance) \n
       \t\t _instance = new N (); \n
       \t\t return _instance; \n
    \t} \n 
 \n
    \t ~N(); \n
    \t void xyz(); \n
 \n
 private: \n 
 static N* _instance;  \n
    \t N();            // no instance of N can be created from the "outside" of N \n
    \t\t\t\t               // protected, if one wants to inherit from that class  \n
    \t N( const N& );              //no new instance via copy constructor  \n
    \t N & operator = (const N &); //no new instance via copy constructor  \n
 };
</CODE>

<strong> Creation of the Data Tables  </strong>

All the following feature are not part of the CRPropa "standard usage". It is merely a set of tools which
were used during the extension of CRPropa to the propagation of nuclei. They haven't been
comment properly and were not tested on different computers with different operational systems
or settings. It's mainly documented for future developer's usage. Redoing the data tables
surely will take some time . . . consider this paragraph as a small head-start.   

As stated above, there are three sets of data files linked with the photo disintegration for:
cross section \f$\sigma\f$ , averaged cross section \f$\sigma_{\rm{avrg}}\f$ or the mean free 
path \f$\lambda\f$. Note, the cross section files should only be used to create the other 
two kind of files, 
never during a simulation. All tools which were used to create these files are contained in
the framework for "advanced usage" and development. It enables the user/coder to e.g. use a different nuclear models 
for the calculation of the photo disintegration cross section for studying systematics. Additionally,
it is possible to decrease the number of tracked channels to gain performance - of course this
is linked with a loss in accuracy. 

Here we use the TALYS framework for the prediction of the photo nuclear cross sections. For that 
a list of isotopes is needed. This can be red from a mysql databases which is defined in 
<CODE> $CRPRopaSrc/MySQLDB/Isotopes/TableDef/CreateIsotopesTable.sql</CODE>. The isotope data
can be filled into the data base by starting 

<CODE>
cd CRPRopaSrcMySQLDB/Isotopes/Fill\n
./FillTable.pl -u 'YourMySQLUser' -p 'YourMySQLPassword' -h 'YourMySQLHost' 
</CODE>


The TALYS calculations were started using another perl script 
<CODE> $CRPropaSrc/GetTalysCrossSections/GetTalysCrossSections.pl </CODE>. After TALYS was properly
installed it was started by calling

<CODE>
cd $CRPropaSrc/GetTalysCrossSections/ \n
./GetTalysCrossSections.pl -EGmmEqualBin "linear" -u 'YourMySQLUser' -p 'YourMySQLPassword' -h 'YourMySQLHost' -ScriptPath '$CRPropaSrc/GetTalysCrossSections' -TALYSDataPath 'YourTALYSDataFolder' -CMin 1 -CMax 26
</CODE>

This starts TALYS for elements with charge numbers CMin<\f$Z\f$<CMax.  The binning in the CMS photon energy
is given to TALYS via the file "gencms" as stored in the same directory. Currently, there are 500 equidistant 
bins from 0.001 MeV  to 250MeV. TALYS was started with the following settings which is given to TALYS via a 
file called "input". 

<CODE>\n
#\n
# General\n
#\n
projectile g\n
element Li\n
mass 7\n
channels y\n
maxchannel 8\n
energy gencms\n
maxN 42\n
maxZ 20\n
</CODE>

These files will be automatically generated. Be aware that the TALYS source file talys.cmb needed to be altered to allow for maxN>32 and maxZ>12.
 
Note, there are two pitfalls if dealing with the TALYS output. TALYS will created very small negative
cross section values - presumably due to numerical inaccuracies. Those values have been set to zero.
Furthermore, there are exclusive channels which are linked with a larger mass- and chargeloss than the mass or the charge
of the nucleus itself. All this channels had a default small and  constant cross section of  \f$10^{-7}\f$mbarn 
These exclusive channels will be sorted out by the GetTalysCrossSections.pl script. Furthermore the channels xs000000.tot
which is not linked with the loss of at least one nucleon are not tracked. For N_11_7 and Li_11_3 TALYS has started
but no exclusive channel files were produced.\n
The current version of TALYS does not start for any helium isotope. It always crashes with the message 
<CODE>TALYS-error: Wrong symbol for element: He</CODE>. This message appears independently of the corresponding 
xml settings: <CODE>element 4</CODE> or <CODE>element He</CODE>. Anyway, the two reasonable stable isotopes 
He-3 and He-4 were independently parametrized similar to what was proposed by Rachen et al.   

<ul>
    <li> <strong> cross section \f$\sigma\f$ files *.cxs: </strong>

    For this one mainly has to extract the cross section information from the TALYS output files
    and create *.cxs files according to the file structure and format as explained above (Note, that the TALYS output
    is not distributed with the CRPropa framework). As parsing through the 
    corresponding files is easier with bash like programming languages, a perl script 
    GetTalysCrossSections/GetTalysCrossSections.pl can perform the *.cxs table creation  
   
    <CODE>
    ./GetTalysCrossSections.pl -CrtTlysDtFlsFtCrPr -TALYSDataPath 'YourTALYSDataFolder'
    </CODE>
   
    The created *.cxs files have to be copied into their expected location at $CRPropaSrc/trunk/TALYSData/. 
    Note, the *.cxs files will be copied into the install directory by the next call of make install.
    
    <CODE>
    cp $CRPropaSrc/GetTalysCrossSections/PD*.cxs $CRPropaSrc/trunk/TALYSData/ \n
    cd $CRPropaSrc/trunk/ \n
    make install
    </CODE>    
    
    <li> <strong> averaged cross section \f$\sigma_{\rm{avrg}}\f$ files *.cax: </strong>

    A simple bash script that starts the creation & thinning (s. below) of these tables can be called
    
    <CODE> 
    cd$CRPropaSrc/examples/CreateReducedTALYSAveragedCrossSectionTable\n 
    ./replaceTALYSData.sh
    </CODE>
    
    The script does not take parameters. If you want to change its settings, you have to modify the script before
    running it. Especially, if you do not intent a table thinning (\f$\alpha=90\%\f$) which by default
    is performed.

    The idea of the "channel thinning" is as follows:\n  
    For each isotope:
    <ol>
    <li>order the i exclusive channels 
    using \f$\sigma_{\rm{avrg}}^{i}\f$  (descending) for each energy,
    <li> calculate the total mean free path \f$\sigma'_{\rm{avrg}}=\sum_i \sigma_{\rm{avrg}}^{i}\f$
         for the isotope, 
    <li> select channels which contribute to the ordered sum 
    \f$ \frac {\sum_i \sigma_{\rm{avrg}}^{n}} {\sigma'_{\rm{avrg}}}>\alpha\f$ in <strong> one </strong> energy bin. 
    Herein, \f$\alpha\f$ is the thinning factor.
    </ol>
    
    In principle, you only have to create the tables once. Then you can created thinned version by using
    the corresponding CRPropa call in the 2nd part of the replaceTALYSData.sh script.
    The script will copy the created or thinned tables into the install directory 
    <CODE> $CRPROPAROOT/share/TabulatedTALYSAveragedCrossSection </CODE>. Thus they will be overwritten
    with the next call of "make install". For permanent usage the tables need to be copied to 
    <CODE> $CRPropaSrc/trunk/TabulatedTALYSAveragedCrossSection/</CODE>. In this case they will be copied
    into the install directory with next call of "make install".
    
    <li> <strong> mean free path \f$\lambda\f$ files *.cmt: </strong>

    Creation and thinning works similar to what is explained for the case of \f$\sigma_{\rm{avrg}}\f$, but here
    one has two scripts - one for the case of the CMB and one for the case of the IRB - in two different directories: 
    <CODE> $CRPropaSrc/examples/CreateReducedTALYSMeanFreePathTable_CMB </CODE> and 
    <CODE> $CRPropaSrc/examples/CreateReducedTALYSMeanFreePathTable_IRB </CODE>. 
    Change to the corresponding directory and simply call.

     <CODE> 
     ./replaceTALYSData.sh
     </CODE>

    One difference to the creation of the  \f$\sigma_{\rm{avrg}}\f$ tables is that in the initiated procedure the 
    minimum and maximum gamma factor of the nuclei as well as the number of bins for the tables have to be given in the script. 
    There is no CRPropa dialog asking for the values as in the case of the \f$\sigma_{\rm{avrg}}\f$ tables.
</ul> 

It is recommended to create a backup of the *.cxs, *.cax and *.cmt files before recreating or changing the table files. 
 
<strong> XML Keywords </strong>

To set up the photo disintegration environment for a simulation, the following
keywords can be used in the xml configuration file for CRPropa.
 
<ul>
 <li> <strong> <NoPhotodisintegration></strong>: by default, the photo disintegration is activated. 
 It can be switched by including this keyword  in the <Interaction> scope. Note, even if photodisintegrations is switched 
 off the mean free path tables will be loaded (but not used).

 <li> <strong> <PhotoDisMFPTabulated></strong>: This tag can be used to select the mean free path \f$\lambda\f$ tables for CMB and IRB. 
 This is the default setting in CRPropa v2.0. 

 <li > <strong> <PhotoDisAvrgdXS></strong>: This keyword can be used to switch CRPropa to use the \f$\sigma_{\rm{avrg}}\f$ tables (undocumented
 feature in this version of CRPropa that might be tested and released with CRPropa 3.0). 
 Using these, the mean free path is calculated by solving the (outer) integral, folding with the photon fields, 
 on the fly. The accuracy of these calculations can be adjusted using the  <NIntegrationCalls> tags. Which defines
 just the number of integrations routine calls and disjoints of the integration range.
 
 It follows the XML-syntax to calculate the mean free path using the IR::IR() IRB background on the fly 
 (Note: IR::IR() returns the Primack-IRB as used in version 1.4 of CRPropa). 

 <CODE>\n
 <Interactions type="Sophia"> \n
  <MaxStep_Mpc value = 10 /> \n
  <NIntegrationCalls value=5 /> \n 
  <PhotoDisAvrgdXS/> \n
 </Interactions> 
 </CODE>
 
 Alternatively the following XML-syntax will calculate the mean free path on the basis of the  variable infrared background TVariableInfrared::TVariableInfrared(...)

 <CODE>\n
 <Interactions type="Sophia"> \n
 <PhotoDisAvrgdXS/> \n 
 <NoIRPhotodisintegration/> \n
 <Photodisintegration_VarIRB/> \n
 <NIntegrationCalls value=5 /> \n
 </Interactions> \n \n
 
 <InfraredBackground type="Variable" > \n
 <SpectralShape type="Shells" /> \n
 <File type="ASCII"> IR_table.data </File> \n
 </InfraredBackground>
 </CODE>

 Here it is important to give the <CODE> <NoIRPhotodisintegration/> </CODE> tag to switch off the usage of IR::IR()!

 \n <STRONG> Important: </STRONG> In version 2.0 of CRPropa the mean free path is hardcoded to be calculated at a redshift of \f$z=0\f$ - even during 
 the on the fly calculations. This is done to perform the mean free path scaling as function of redshift as shortly described in TIRBzEvolutionModel 
 and in the CRPropa version 2.0 paper. Thus, a variable infrared model with a redshift dependency would always be evaluated at $z=0$ while calculating 
 the mean free path.    
 </ul>
 
 <strong> Todos & Notes </strong>
 
 <ul>
 <li> The reading implementation of the data files (*.cxs, *.cax, *.cmt) is similar and can be done in one single routine in TabulatedTALYSY which is then inherited to the concrete classes (e.g. TabulatedTALYSMeanFreePath). 
A similar reduction of code should be possible with e.g. the CreateTables(. . . ) and ReduceTables(. . . ) functions.
 <li> The names of variables and functions often include the term \"energy\" which is somewhat misleading as one deals with a tabulated version of a general function \f$Y(x)\f$. Here \"energy\" simply corresponds to the \f$x\f$ value of the function. A consistent renaming is recommended to make the code more readable. 
 <li> Further checks are needed if there are more stable light nuclei which have to be included.
 <li> Allow for the combined usage of mean free path tables and on the fly calculation using the averaged cross
 section. E.g. CMB from mean free path table and variable IR from on the fly calculations.
 <ul/>
 
*/




using namespace std;

class TabulatedTALYSY {
 public:

  ~TabulatedTALYSY();
 
  void ASCIIOutExclY(unsigned long int InputNucleusId,
			     unsigned long int ExclusiveChannelId);
  /**< Creates an ASCII file which holds x, Y(x) as stored in the table for the "exclusive" channel: ExclusiveChannelId
     for nucleus InputNucleusId.*/

  void RootPlotExclY(unsigned long int InputNucleusId,
			     unsigned long int ExclusiveChannelId);
  /**< Creates a root plot which holds x, Y(x) as stored in the table for a "exclusive" channel:  ExclusiveChannelId
     for nucleus InputNucleusId.*/
  
  inline double TotalYValueBelowFXMin(unsigned long int InputNucleusId,
				       unsigned long int OutputNucleusId){return(0.);}
  /**< Defines which Y(x) value is returned if x<FXmin for total channel: InputNucleusId-> OutPutNucleusId (out of table range).*/
  
  inline double TotalYValueAboveFXMax(unsigned long int InputNucleusId,
				       unsigned long int OutputNucleusId){return(0.);}
  /**< Defines which Y(x) value is returned if x>FXmax for total channel: InputNucleusId-> OutPutNucleusId (out of table range).*/
  
  inline double YValueNotAvailable(){return(0.);}
  /**< Defines which Y(x) value is returned if no entry is available.*/

  inline double ExclYValueBelowFXMin(unsigned long int InputNucleusId,
				      unsigned long int ExclusiveChannelId){return(0.);}
  /**< Defines which Y value is returned if x<FXmin for exclusive channel: ExclusiveChannelId for nucleus InputNucleusId (out of table range).*/
  
  virtual double ExclYValueAboveFXMax(unsigned long int InputNucleusId,
				      unsigned long int ExclusiveChannelId){return(0.);}
  /**< Defines which Y value is returned if x>FXmax for exclusive channel: ExclusiveChannelId for nucleus InputNucleusId (out of table range).*/
  
  inline double ExclSumYValueBelowFXMin(unsigned long int InputNucleusId){return(0.);}
  /**< Defines which "summed" Y value is returned if x<FXmin for exclusive channel:  nucleus InputNucleusId (out of table range). "Summed" means the sum over all exclusive channels.*/
  
  virtual double ExclSumYValueAboveFXMax(unsigned long int InputNucleusId){return(0.);}
  /**< Defines which "summed" Y value is returned if x>FXmax for exclusive channel:  nucleus InputNucleusId (out of table range). "Summed" means the sum over all exclusive channels.*/
  
  virtual double GetExclusiveChannelY(unsigned long int InputNucleusId,
				      int ExclusiveChannelId,
				      double Energy);
  /**< Returns the Y(x) for the exclusive channel ExclusiveChannelId for InputNucleusId using linear approximation between the tabulated data.*/
  
  virtual double GetExclusiveSummedY(unsigned long int InputNucleusId,
				     double Energy);
  /**< Returns the "summed" Y(x) for InputNucleusId using linear approximation between the tabulated data. "Summed" means the sum over all exclusive channels.*/
  
  virtual int GetEnergyBinIdBelow(double Energy);
  /**< Returns the Id of the closest bin which is closest but below Energy. */
  
  virtual int GetEnergyBinIdAbove(double Energy);
  /**< Returns the Id of the closest bin which is closest but above Energy. */
  
  virtual double GetTabExclChannelYValue(unsigned long int InputNucleusId,
					 unsigned long int ExclusiveChannelId,
					 int EnergyBin);
  /**< Get the tabulated Y value for the energy bin EnergyBin for the exclusive channel ExclusiveChannelId for nucleus  InputNucleusId (start to count bin at 0!). */
  
  virtual double GetTabExclSummedYValue(unsigned long int InputNucleusId,
					int EnergyBin);
  /**<  Get the "summed" tabulated Y values for the energy bin EnergyBin for the exclusive channel ExclusiveChannelId for nucleus  InputNucleusId (start to count bin at 0!). "Summed" means the sum over all exclusive channels.*/
  
  virtual double GetEnergyForEnergyBinId(unsigned long int EnergyBin);
  /**< Get the energy for energy bin EnergyBin bin.*/
  
  std::vector<unsigned long int> GetListOfInputNuclei();
  /**< Returns a vector of valid input nuclei.*/
  
  virtual vector<unsigned long int> GetListOfAllExclusiveChannels(unsigned long int InputNucleusId); 			    
  /**< Returns a vector of possible exclusive channels for nucleus InputNucleusId.*/


  inline double GetMaximumXValue(){ return (fXMax);} ;
  /**< The tabulated values are given for [fXMin,fXMax]. This returns fXMax */	 
  
  
  inline double GetMinimumXValue(){ /*std::cout<<"GetMinimum()"<<std::endl;*/ return (fXMin);} ;
  /**< the tabulated values are given for [fXMin,fXMax]. This returns fXMax */	 
  
  inline int GetNXBins(){ return (fXNSteps);} ;
  /**< return the number of x bins in the tables. */	 
  
  std::vector<int> GetDMassAndDChargeExclId  (int ExclusiveChannelId);
  /**<Calculate the difference in mass and charge number which is linked with exclusive channel ExclusiveChannelId:
     Returns a vector: 1. component deltaMass 2. component deltaCharge*/

  std::vector<unsigned long int> GetMassAndChargeNuclId  (unsigned long int NucleusId);
  /**<Converts NucleiId back into mass- A and chargenumber Z. 
     Returns a vector: 1. component MassNumber 2. component ChargeNumber*/
  
  std::vector<int> GetDMassAndDChargeNuclId(unsigned long int NucleusId1,
					    unsigned long int NucleusId2);
  /**<Calculate the difference in mass and charge number which is linked with 
     nucleus transition:  NucleusId1 ->  NucleusId2.
     Returns a vector: 1. component deltaMass 2. component deltaCharge*/

  virtual void RootFileWithAllReactions(double fXMin,
					double fXMax,
					double fXNSteps,
					double redshift);
  /**<Create a set of root files which contain plots of the data as stored in the tables. 
     It uses the three functions that follow.*/
  
  virtual TGraph* RootPlotYVsX(unsigned long int InputNucleusId,
			       long unsigned int OutId,
			       double fXMin,
			       double fXMax,
			       double fXNSteps,
			       double redshift,
			       int TotalSumOrExclusive);
  /**< Plot routine.*/
  
  virtual TGraph* RootPlotExclSumYVsX(unsigned long int InputNucleusId,
				      double fXMin,
				      double fXMax,
				      double fXNSteps,
				      double redshift);
  /**< Plot routine.*/

  virtual void CreateTables(std::vector<PhotonBackground*> PhotonBackgrounds_Vec,
			    std::vector<TabulatedTALYSY*> tTabulatedTALYSY_Vec,
			    PD_TabDataType tPD_TabDataType,
			    TVector3D Position,
			    double epsilon0, 
			    double epsilon_min, 
			    double epsilon_max,
			    std::string folder,
			    int NIntegratorCallsdouble,
			    double EMinMFPTableVal,
			    double EMaxMFPTableVal,
			    int NBinsMFPTableVal){throw TCrpErr("TabulatedTALYSY::CreateTables(1) called but not implemented");}
  
  virtual void CreateTables(std::string folder){throw TCrpErr("TabulatedTALYSY::CreateTables(2) called but not implemented");}

  
  /**< Returns an value of type enum PD_TabDataType {CMB_MFPTab, IRB_MFPTab} which defines if the tabulated MFP is linked with CMB or IRB.*/
  PD_TabDataType GetPD_TabDataType(){return fPD_TabDataType;}
  
  std::vector<double> fExclusiveYData_Vec;
  std::vector<double> fExclusiveSumYData_Vec;
  std::vector< std::vector< unsigned long int >  > fTALYSExclYId_Matrix;
  std::vector< std::vector< unsigned long int >  > fTALYSInitNucleusId_Matrix;
  
  //define the ranges and stepsize of the TALYS data in the folder TALYSData
  double fXMin;        //in MeV allowed range is (fMinPhotonEn,fMaxPhotonEn]
  double fXMax;        //in MeV
  int fXNSteps;
  double fXStepsize;   //in MeV
  double fthinning_alpha;

  double fPDFileType;
  std::string fXUnit; 
  std::string fXDescription; 
  std::string fYUnit; 
  std::string fYDescription; 

 protected:
  TabulatedTALYSY();  
  PD_TabDataType fPD_TabDataType;

 private:
  TabulatedTALYSY(const TabulatedTALYSY&);
  TabulatedTALYSY & operator = (const TabulatedTALYSY&);
};

#endif

