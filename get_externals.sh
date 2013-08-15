#!/bin/sh

CURDIR=`pwd`

# CFITSIO_VERSION=3006
CFITSIO_VERSION=3300
#CLHEP_VERSION=2.0.1.1
CLHEP_VERSION=2.0.4.3

#CFITSIO_SRC=http://apcauger.in2p3.fr/CRPropa/cfitsio$CFITSIO_VERSION.tar.gz
#CLHEP_SRC=http://apcauger.in2p3.fr/CRPropa/clhep-$CLHEP_VERSION.tgz

# CFITSIO_SRC=http://astro.uni-wuppertal.de/~nils/cfitsio$CFITSIO_VERSION.tar
CFITSIO_SRC=ftp://heasarc.gsfc.nasa.gov/software/fitsio/c/cfitsio$CFITSIO_VERSION.tar.gz
#CLHEP_SRC=http://astro.uni-wuppertal.de/~nils/clhep-$CLHEP_VERSION.tar
CLHEP_SRC=http://proj-clhep.web.cern.ch/proj-clhep/DISTRIBUTION/tarFiles/clhep-$CLHEP_VERSION.tgz

WGET_OK=0
 
test -n "`which curl 2> /dev/null`" && WGET_OK=1 && WGET="curl -O"
test -n "`which wget 2> /dev/null`" && WGET_OK=1 && WGET="wget"

cd External

if test $WGET_OK -eq 1 
then
  $WGET $CLHEP_SRC
  $WGET $CFITSIO_SRC
#  echo
else
  echo "Nor wget neither curl was found on this system. Cannot get automatically external program tar balls. Sorry, exiting." >&2 && exit 1;
fi

GCC_OK=0
test -n "`which g++ 2> /dev/null`" && GCC_OK=1

if test $GCC_OK -ne 1
then
  echo "No Gnu compiler (gcc) found on the system. Exiting." >&2 && exit 2;
fi

GCC_VERSION=`gcc -dumpversion | awk -F"." '{print $1}'`
if test $GCC_VERSION -lt 3
then
  echo "GCC version is lower than 3. Cannot continue. Exiting." >&2 && exit 3;
fi

CC=gcc
export CC
CXX=g++
export CXX

gzip -d `basename $CFITSIO_SRC`
tar xvf `basename $CFITSIO_SRC | sed 's/\.gz//' ` 
gzip -d `basename $CLHEP_SRC`
tar xvf `basename $CLHEP_SRC | sed 's/\.tgz/\.tar/' `

# -- how many procs
NB_PROC=1
if test `uname` = "Linux"
then
  NB_PROC=`grep processor < /proc/cpuinfo | wc -l | awk '{print $1}'`
fi
if test `uname` = "Darwin"
then
  NB_PROC=`sysctl hw.ncpu | awk '{print $2}'`
fi
test $NB_PROC -ge 2 && MAKE_OPT="-j $NB_PROC"

# -- cfits io compil
cd cfitsio
./configure --prefix=`pwd` && make $MAKE_OPT && make install
cd $CURDIR/External 

cd $CLHEP_VERSION/CLHEP
./configure --prefix=$CURDIR/External/ && make $MAKE_OPT && make install

cd $CURDIR

#echo 
#echo "####################################################"
#echo "Now you should set CLHEP_DIR=$CURDIR/External/ and "
#echo "CFITSIO_DIR=$CURDIR/External/cfitsio to continue the "
#echo "installation of CRPropa..."
#echo "####################################################"
#echo 

