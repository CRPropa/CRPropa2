# Simple to script to download and copy the cross section data for photo disintegration
# author: Nils Nierstenhoefer (nierstenhoefer@physik.uni-wuppertal.de)
#!/bin/sh

CURDIR=`pwd`

echo -e "\n\n\n**********************************************************************************************\n"
echo -e "*                                                                                            *\n"
echo -e "*              INSTALLATION of CRPropa Photo Disintegration Data Package                     *\n"
echo -e "*                                                                                            *\n"
echo -e "**********************************************************************************************\n"

if [ -e "CrossSectionPackage2_6_Kneiske.tar.gz" ]; then
	echo -e "\nFile CrossSectionPackage2_6_Kneiske.tar.gz exists. Replace it."
	rm CrossSectionPackage2_6_Kneiske.tar.gz
else 
	echo -e "\n Preparing for download."
fi 

PDCrossSec_SRC=https://crpropa.desy.de/images/2/2f/CrossSectionPackage2_6_Kneiske.tar.gz
#PDCrossSec_SRC=http://astro.uni-wuppertal.de/~nils/CrossSectionPackage2_6_Kneiske.tar.gz

test -n "`which curl 2> /dev/null`" && WGET_OK=1 && WGET="curl -O"
test -n "`which wget 2> /dev/null`" && WGET_OK=1 && WGET="wget --no-check-certificate"

echo -e "\n1. Try to download photo disintegration cross section data package . . . \n"

if test $WGET_OK -eq 1 
then
  $WGET $PDCrossSec_SRC
  echo ""
else
  echo -e "Nor wget neither curl was found on this system. Cannot get automatically external program tar balls. Sorry, exiting." >&2 && exit 1;
fi


echo -e "\n2. Try to unpack photo disintegration cross section package . . . "

tar -xzvf CrossSectionPackage2_6_Kneiske.tar.gz

echo -e "\n3. Try to create directory structure . . . "
if [ -d "TabulatedTALYSMeanFreePath" ]; then
	echo -e "\n Directory TabulatedTALYSMeanFreePath exists"
else 
	echo -e "\n Directory TabulatedTALYSMeanFreePath does not exists"
	mkdir TabulatedTALYSMeanFreePath
fi 

if [ -d "TabulatedTALYSMeanFreePath/CMB" ]; then
	echo -e "\n Directory TabulatedTALYSMeanFreePath/CMB exists"
else 
	echo -e "\n Directory TabulatedTALYSMeanFreePath/CMB does not exists"
	mkdir TabulatedTALYSMeanFreePath/CMB
fi 

if [ -d "TabulatedTALYSMeanFreePath/IRB" ]; then
	echo -e "\n Directory TabulatedTALYSMeanFreePath/IRB exists"
else 
	echo -e "\n Directory TabulatedTALYSMeanFreePath/IRB does not exists"
	mkdir TabulatedTALYSMeanFreePath/IRB
fi 

#if [ -d "TabulatedTALYSMeanFreePath/CMB_IRB" ]; then
#	echo -e "\n Directory TabulatedTALYSMeanFreePath/CMB_IRB exists"
#else 
#	echo -e "\n Directory TabulatedTALYSMeanFreePath/CMB_IRB does not exists"
#	mkdir TabulatedTALYSMeanFreePath/CMB_IRB
#fi 

if [ -d "TALYSData" ]; then
	echo -e "\n Directory TALYSData exists"
else 
	echo -e "\n Directory TALYSData does not exists"
	mkdir TALYSData
fi 

if [ -d "TabulatedTALYSAveragedCrossSection" ]; then
	echo -e "\n Directory TabulatedTALYSAveragedCrossSection exists"
else 
	echo -e "\n Directory TabulatedTALYSAveragedCrossSection does not exists"
	mkdir TabulatedTALYSAveragedCrossSection
fi 

echo -e "\n3. Copying TALYS Data into source folder . . . "
cp CrossSectionPackage2_6_Kneiske/TALYSData/PD*.cxs TALYSData/

echo -e "\n4. Copying TALYS averaged cross section data into source folder . . . "
cp CrossSectionPackage2_6_Kneiske/TabulatedTALYSAveragedCrossSection_XSBins500_ABins200_thn90/PD*.cax TabulatedTALYSAveragedCrossSection/

echo -e "\n5. Copying TALYS mean free path data into source folder . . . "
cp -r CrossSectionPackage2_6_Kneiske/TabulatedTALYSMeanFreePath_XSBins500_ABins200_MFPBins200_thn90/* TabulatedTALYSMeanFreePath/

echo -e "\n\t . . . data was copied to source directory. They will only ready for usage   A F T E R   the next call of \"make install\"!!!" 

echo -e "\n\n\nNote, the recommended tables with a thinning factor of 90% have been installed. "
echo -e "In case you want to use a different thinning check the doxygen page for further documentation.\n\n\n"

