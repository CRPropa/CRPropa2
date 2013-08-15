SUB_CPPFLAGS = \
	-I@abs_top_builddir@/sophia/src \
	-I@abs_top_builddir@/sibyll \
	-I@abs_top_builddir@/dint/src   \
	-I@abs_top_builddir@/External/tinyxml \
	-I@abs_top_builddir@/src \
	-I@abs_top_builddir@/src/Algo \
	-I@abs_top_builddir@/src/Config \
	-I@abs_top_builddir@/src/Exceptions \
	-I@abs_top_builddir@/src/Interactions \
	-I@abs_top_builddir@/src/Particles \
	-I@abs_top_builddir@/src/Environment \
	-I@abs_top_builddir@/src/Gas \
	-I@abs_top_builddir@/src/MagField \
	-I@abs_top_builddir@/src/Observers \
	-I@abs_top_builddir@/src/Sources \
	-I@abs_top_builddir@/src/PhotonBackground \
	-I@abs_top_builddir@/sophia \
	-DDEFAULT_IRBzRedshiftEvol_Kneiske=\"@datadir@/IRBzRedshiftEvol_Kneiske/\" \
	-DDEFAULT_F77_PROTON_DIR=\"@datadir@/proton_f77/\" \
	-DDEFAULT_SOPHIA_DIR=\"@datadir@/proton_sophia/\"  \
	-DDEFAULT_SHOWER_DIR=\"@datadir@/dint_tables/\"    \
	-DDEFAULT_HLT_DIR=\"@datadir@/HalfLifeTable/\"     \
	-DDEFAULT_TALYS_DIR=\"@datadir@/TALYSData/\"	\
	-DDEFAULT_TALYSAveraged_DIR=\"@datadir@/TabulatedTALYSAveragedCrossSection/\"	\
	-DDEFAULT_TALYSTabMeanFreePath_DIR=\"@datadir@/TabulatedTALYSMeanFreePath/\"	
#	-I@abs_top_builddir@/src/IRBackground 


SUB_LNKFLAGS = \
	-LConfig -lConfig \
	-LInteractions -lInteractions \
	-LParticles -lParticles \
	-LEnvironment -lEnvironment \
	-LGas -lGas \
	-LMagField -lMagField \
	-LObservers -lObservers \
	-LSources -lSources \
	-LPhotonBackground -lPhotonBackground \
	-L@abs_top_builddir@/sophia -lsophia \
	@abs_top_builddir@/sibyll/libsibyll.a \
	-L@abs_top_builddir@/dint -ldint \
	-L@abs_top_builddir@/External/tinyxml -ltinyxml\
	-lfftw3
#	-LIRBackground -lIRBackground 

export SUB_CPPFLAGS
export SUB_LNKFLAGS 
export FFLAGS

AM_CPPFLAGS = $(GNL_CPPFLAGS)

#SUBDIRS = External doc sophia dint src 
SUBDIRS = External sophia sibyll dint src  doc

install-data-local:
	test -d @datadir@ || mkdir @datadir@
	test -d @datadir@/tables || mkdir -p @datadir@/tables
	cp -rpu `ls @abs_top_builddir@/dint/src/tables/* | grep -v -E .svn` @datadir@/tables/.
	test -d @datadir@/dint_tables || mkdir -p @datadir@/dint_tables
	mv -u @datadir@/tables/* @datadir@/dint_tables/
	test -d @datadir@/proton_f77 || mkdir -p @datadir@/proton_f77
	cp -rpu `ls @abs_top_builddir@/src/Interactions/proton_f77/* | grep -v -E .svn`  @datadir@/proton_f77/.
	test -d @datadir@/proton_sophia || mkdir -p @datadir@/proton_sophia
	cp -rpu `ls @abs_top_builddir@/src/Interactions/proton_sophia/* | grep -v -E .svn`  @datadir@/proton_sophia/.
#	test -d @datadir@/Crptests || mkdir -p @datadir@/Crptests
#	cp -rpu `ls @abs_top_builddir@/Crptests/* | grep -v -E .svn` @datadir@/Crptests/.
	test -d @datadir@/HalfLifeTable || mkdir -p @datadir@/HalfLifeTable
	cp -rpu `ls @abs_top_builddir@/HalfLifeTable/* | grep -v -E .svn` @datadir@/HalfLifeTable/.		
	cp -rpu `ls @abs_top_builddir@/HalfLifeTable/* | grep -v -E .svn` @datadir@/HalfLifeTable/.			
	test -d @datadir@/IRBzRedshiftEvol_Kneiske || mkdir -p @datadir@/IRBzRedshiftEvol_Kneiske
	cp -rpu `ls @abs_top_builddir@/IRBzRedshiftEvol_Kneiske/* | grep -v -E .svn` @datadir@/IRBzRedshiftEvol_Kneiske/.		
	cp -rpu `ls @abs_top_builddir@/IRBzRedshiftEvol_Kneiske/* | grep -v -E .svn` @datadir@/IRBzRedshiftEvol_Kneiske/.			

#	cp -rpu @abs_top_builddir@/TALYSData @datadir@/
#	cp -rpu @abs_top_builddir@/TabulatedTALYSAveragedCrossSection @datadir@/	
#Copy TALYS data into install directory 
	test -d @datadir@/TALYSData || mkdir -p @datadir@/TALYSData	
	(test -d  @abs_top_builddir@/TALYSData/ ) || touch  @abs_top_builddir@/TALYSData/buffer
	cp -rpu `ls @abs_top_builddir@/TALYSData/* | grep -v -E .svn` @datadir@/TALYSData/.		
#Copy TALYS averaged data into install directory 
	test -d @datadir@/TabulatedTALYSAveragedCrossSection || mkdir -p @datadir@/TabulatedTALYSAveragedCrossSection	
	(test -d  @abs_top_builddir@/TabulatedTALYSAveragedCrossSection/ ) || touch  @abs_top_builddir@/TabulatedTALYSAveragedCrossSection/buffer
	cp -rpu `ls @abs_top_builddir@/TabulatedTALYSAveragedCrossSection/* | grep -v -E .svn` @datadir@/TabulatedTALYSAveragedCrossSection/.		
#Copy TALYS mean free path data into install directory 
	test -d @datadir@/TabulatedTALYSMeanFreePath || mkdir -p @datadir@/TabulatedTALYSMeanFreePath
	(test -d  @abs_top_builddir@/TabulatedTALYSMeanFreePath/) || touch @abs_top_builddir@/TabulatedTALYSMeanFreePath//buffer

#CMB_IRB
#	test -d @datadir@/TabulatedTALYSMeanFreePath/CMB_IRB || mkdir -p @datadir@/TabulatedTALYSMeanFreePath/CMB_IRB
#	(test -e  @abs_top_builddir@/TabulatedTALYSMeanFreePath/CMB_IRB/PD*.cmt) || touch @abs_top_builddir@/TabulatedTALYSMeanFreePath/CMB_IRB/buffer
#
#	cp -rpu `ls @abs_top_builddir@/TabulatedTALYSMeanFreePath/CMB_IRB/* | grep -v -E .svn` @datadir@/TabulatedTALYSMeanFreePath/CMB_IRB/.	


#IRB
	test -d @datadir@/TabulatedTALYSMeanFreePath/IRB || mkdir -p @datadir@/TabulatedTALYSMeanFreePath/IRB
	(test -d  @abs_top_builddir@/TabulatedTALYSMeanFreePath/IRB/) || touch @abs_top_builddir@/TabulatedTALYSMeanFreePath/IRB/buffer

	cp -rpu `ls @abs_top_builddir@/TabulatedTALYSMeanFreePath/IRB/* | grep -v -E .svn` @datadir@/TabulatedTALYSMeanFreePath/IRB/.	

#CMB
	test -d @datadir@/TabulatedTALYSMeanFreePath/CMB || mkdir -p @datadir@/TabulatedTALYSMeanFreePath/CMB
	(test -d  @abs_top_builddir@/TabulatedTALYSMeanFreePath/CMB/) || touch @abs_top_builddir@/TabulatedTALYSMeanFreePath/CMB/buffer
#	(test -e  @abs_top_builddir@/TabulatedTALYSMeanFreePath/CMB/PD*.cmt) || touch @abs_top_builddir@/TabulatedTALYSMeanFreePath/CMB/buffer

	cp -rpu `ls @abs_top_builddir@/TabulatedTALYSMeanFreePath/CMB/* | grep -v -E .svn` @datadir@/TabulatedTALYSMeanFreePath/CMB/.		

EXTRA_DIST = sophia/src/sophia.h sibyll/sibyll.h src/Interactions/proton_f77 src/Interactions/proton_sophia dint/src/tables sysdep.h get_externals.sh GetPDCrossSections.sh HalfLifeTable doc/Doxyfile COPYRIGHT doc/UserGuide.tex doc/precision_traj.ps doc/PhotonDensityCollection2.eps doc/crpdetector.eps doc/AllIRBzEvolutionModelsCan.eps IRBzRedshiftEvol_Kneiske examples/GettingStarted Plot

dist-hook:
	rm -rf `find $(distdir) -type d -name CVS`
	rm -rf `find $(distdir) -type d -name ".svn"`
