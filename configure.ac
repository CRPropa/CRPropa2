#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
m4_define(crpropa_version_major, 2)
m4_define(crpropa_version_minor, 0)

AC_PREREQ(2.59)
AC_INIT([CRPropa],[crpropa_version_major().crpropa_version_minor()],[crpropa_support@apc.univ-paris7.fr])
AM_CONFIG_HEADER([config.h])
AM_INIT_AUTOMAKE
LT_INIT
AC_SUBST([LIBTOOL_DEPS])
AC_CONFIG_MACRO_DIR([m4])
ACLOCAL_AMFLAGS='-I m4'

CRPROPA_VERSION_MAJOR=crpropa_version_major()
CRPROPA_VERSION_MINOR=crpropa_version_minor()

AC_SUBST(CRPROPA_VERSION_MAJOR)
AC_SUBST(CRPROPA_VERSION_MINOR)

AC_COPYRIGHT([
	Copyright 2010.
	All Rights Reserved.
	This file is part of the open source package CRPropa which is distributed
	under the Gnu LGPL.])

# Where should everything be installed by default?  Here, we want it
# to be installed directly in 'bin', 'lib', 'include' subdirectories
# of the directory where configure is run.  The default would be
# /usr/local.
AC_PREFIX_DEFAULT([`pwd`])

# Checks for programs.
AC_LANG([C++])
AC_PROG_CXX([g++ clang++])
AC_PROG_CC([gcc clang])
AC_PROG_F77([gfortran ifort g77 g95 f77])
AC_PROG_INSTALL([/usr/bin/install])

# Use libraries in tools
AC_PROG_LIBTOOL
AC_SUBST(LIBTOOL_DEPS)
AC_LTDL_DLLIB

# Add mixing f77 with C feature
AC_F77_LIBRARY_LDFLAGS

# Checks for external libraries.
LDFLAGS_INI="$LDFLAGS"

# -- CLHEP Test --
  AC_ARG_WITH(
	  [clhep-path],
	  AS_HELP_STRING(
	    [--with-clhep-path=DIR],
			[clhep-config binary is in DIR]
		)
	)

	AC_PATH_PROG(
	  [CLHEP_CONFIG],
	  [clhep-config],
		[notfound],
		[$with_clhep_path:$PWD/External/bin:$PATH])
	if test "$CLHEP_CONFIG" == "notfound"
	then
	  AC_MSG_ERROR([clhep-config not found. See ./configure --help to know how to set a proper clhep path])
  fi
	# Test version clhep. >= 2.0.0.0 (tested version 2.0.1.1)
  CLHEP_CONFIG_VERSION=`$CLHEP_CONFIG --version | cut -f2 -d' '`
  if test "`echo $CLHEP_CONFIG_VERSION | cut -f1 -d.`" != "2"
  then 
	    AC_MSG_ERROR([Need clhep version 2. Found $CLHEP_CONFIG_VERSION.])
  else 
    CPPFLAGS="$CPPFLAGS `$CLHEP_CONFIG --include`"
    LDFLAGS_CLHEP="$LDFLAGS_INI `$CLHEP_CONFIG --libs`"
  fi
  LDFLAGS="$LDFLAGS_CLHEP"

  AC_CHECK_HEADERS(
    [CLHEP/Random/RandGauss.h CLHEP/Units/SystemOfUnits.h CLHEP/Vector/ThreeVector.h],
		[AC_MSG_RESULT([check for clhep include files... yes])],
		[AC_MSG_ERROR([Cannot find some needed CLHEP include files])]
	)

	AC_LINK_IFELSE(
	  [AC_LANG_SOURCE([[
	#include "CLHEP/Random/RandGauss.h"
	using namespace CLHEP;
	int main () { 
	  RandGauss::shoot() ; 
		return 0;
	}
	]]
	)],
	  [AC_MSG_RESULT([check for clhep link... yes])],
		[AC_MSG_ERROR([Cannot link with CLHEP lib. Check CLHEP installation.])]
	)

# -- ROOT Test --
  AC_ARG_WITH(
	  [root],
		AS_HELP_STRING(
		  [--with-root(=DIR)],
			[using Root output option, assuming Root installed in DIR]
		)
	)

  AM_CONDITIONAL(ROOTCOMPIL,test "x$with_root" != "x")

#  if test "x$with_root" != "x" 
#	then
 	  # Using root
    AC_PATH_PROG(
	    [ROOT_CONFIG],
  		[root-config],
  		[notfound],
  		[$with_root/bin:$ROOTSYS/bin:$PATH]
  	)
  	if test "$ROOT_CONFIG" == "notfound"
  	then
  	  AC_MSG_ERROR([root-config not found. See ./configure --help to know how to set a proper root path])
  	else
      CPPFLAGS="$CPPFLAGS `$ROOT_CONFIG --cflags`"
			LDFLAGS_ROOT=" `$ROOT_CONFIG --libs` "
	  fi
#	fi
	
# -- CFITSIO Test --
  AC_ARG_WITH(
	  [cfitsio-include],
    AS_HELP_STRING(
			  [--with-cfitsio-include=DIR],
				[cfitsio library include files are in DIR]
		)
	)
  if test "x$with_cfitsio_include" != "x" 
	then
    CPPFLAGS="$CPPFLAGS -I$with_cfitsio_include"
  else
	  CPPFLAGS="$CPPFLAGS -I$PWD/External/cfitsio/"
	fi
	
# Test include fitsio.h found in CPPFLAGS
#AC_PREPROC_IFELSE(
#  [AC_LANG_PROGRAM([[const float v = CFITSIO_VERSION ;]])],
#	[AC_MSG_RESULT([checking for fitsio.h... yes])],
#	[AC_MSG_ERROR([Cannot find fitsio.h. See ./configure --help to know how to set a proper cfitsio include directory])] 
#)
AC_CHECK_HEADER(
  [fitsio.h],
	[],
	[AC_MSG_ERROR([Cannot find fitsio.h. See ./configure --help to know how to set a proper cfitsio include directory])]
)

LDFLAGS="$LDFLAGS_INI"
# Test CFITSIO_VERSION >= 3.006
  AC_RUN_IFELSE(
	  AC_LANG_SOURCE([[
#include <iostream>
#include "fitsio.h"
int main() { 
  int r;
	r = !( CFITSIO_VERSION >= 3. ) ;
	if (r) std::cerr << "cfitsio version = " << CFITSIO_VERSION << std::endl; 
	return r ;
}
	  ]]),
	  [AC_MSG_RESULT([checking for cfistio version... yes])],
	  [AC_MSG_ERROR([Bad cfitsio version. Need >= 3.])]
	)
	
	AC_ARG_WITH(
	  [cfitsio-library],
 		AS_HELP_STRING(
		  [--with-cfitsio-library=DIR],
			[cfitsio library file is in DIR]
		)
	)
	if test "x$with_cfitsio_library" != "x"
	then
	  LDFLAGS_CFITSIO="$LDFLAGS -L$with_cfitsio_library"
	else
	  LDFLAGS_CFITSIO="$LDFLAGS -L$PWD/External/cfitsio/lib"
	fi

  LDFLAGS="$LDFLAGS_CFITSIO"
	AC_CHECK_LIB(
	  [cfitsio],
		[ffvers],
		[],
		[AC_MSG_ERROR(
		  [Cannot find cfitsio library. See ./configure --help to know how to set a proper library path]
		)]
  )

# -- fftw3 Test --
  AC_ARG_WITH(
	  [fftw-include],
    AS_HELP_STRING(
			  [--with-fftw-include=DIR],
				[FFTW3 library include files are in DIR]
		)
	)
  if test "x$with_fftw_include" != "x" 
	then
    CPPFLAGS="$CPPFLAGS -I$with_fftw_include"
	fi

  if test "x$with_fftw_include" == "x"	
     then 
     PKG_CHECK_MODULES([fftw3], [fftw3])
     fi   	
AC_CHECK_HEADER(
  [fftw3.h],
	[],
	[AC_MSG_ERROR([Cannot find fftw3.h. See ./configure --help to know how to set a proper fftw3 include directory])]
)

LDFLAGS="$LDFLAGS_INI"
	
	AC_ARG_WITH(
	  [fftw-library],
 		AS_HELP_STRING(
		  [--with-fftw-library=DIR],
			[FFTW3 library file is in DIR]
		)
	)
	if test "x$with_fftw_library" != "x"
	then
	  LDFLAGS_FFTW="$LDFLAGS -L$with_fftw_library"
	fi

  LDFLAGS="$LDFLAGS_INI $LDFLAGS_CLHEP $LDFLAGS_CFITSIO $LDFLAGS_FFTW $LDFLAGS_ROOT $FLIBS"

# Checks for header files.
	AC_HEADER_STDC
	AC_CHECK_HEADERS([arpa/inet.h malloc.h netinet/in.h stdlib.h string.h unistd.h assert.h ctype.h fftw3.h])
  AC_CHECK_HEADERS([iostream fstream sstream string vector queue cfloat functional])
	AC_CHECK_HEADERS([TROOT.h])

# Checks for typedefs, structures, and compiler characteristics.
	AC_HEADER_STDBOOL
	AC_C_CONST
	AC_C_INLINE
	AC_TYPE_SIZE_T

# Checks for library functions.
#	AC_FUNC_MALLOC
	AC_FUNC_STAT
AC_CHECK_FUNCS([floor pow sqrt strchr strtoul])

# Checks for programs to build docs
AC_CHECK_PROG(LATEX, latex, latex)
AC_CHECK_PROG(DVIPS, dvips, dvips)
AC_CHECK_PROG(DVIPDF, dvipdf, dvipdf)
AC_CHECK_PROG(DOXYGEN, doxygen, doxygen)

AM_CONDITIONAL(LATEXTEST,test x$LATEX = xlatex)
AM_CONDITIONAL(DVIPSTEST,test x$DVIPS = xdvips)
AM_CONDITIONAL(DVIPDFTEST,test x$DVIPDF = xdvipdf)
AM_CONDITIONAL(DOXYGENTEST,test x$DOXYGEN = xdoxygen)

# global CPP flags
AC_SUBST(CPPFLAGS)

# global CXX flags
AC_SUBST(CXXFLAGS)

FFLAGS_INI="$FFLAGS"
#if test "x$withpic" != "x"
#   then
FFLAGS="-fPIC $FFLAGS_INI"
#fi

AC_SUBST(FFLAGS)

# Make makefiles
#			src/IRBackground/Makefile
	AC_CONFIG_FILES([External/Makefile
			External/tinyxml/Makefile
			Makefile
			dint/Makefile
			sophia/Makefile
			sibyll/Makefile
			src/Makefile
			src/Config/Makefile
			src/Environment/Makefile
			src/Interactions/Makefile
			src/Gas/Makefile
			src/MagField/Makefile
			src/Observers/Makefile
			src/Particles/Makefile
			src/Sources/Makefile
			src/PhotonBackground/Makefile
			doc/Makefile
			])
	AC_OUTPUT

