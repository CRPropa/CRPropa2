
AM_CPPFLAGS = $(SUB_CPPFLAGS)
LDADD  = $(SUB_LNKFLAGS) 

export SUB_CPPFLAGS

SUBDIRS = Config Environment Gas Interactions MagField Observers Particles Sources PhotonBackground

bin_PROGRAMS = CRPropa
CRPropa_LDADD = $(SUB_LNKFLAGS) 

CRPropa_SOURCES = CRPropa.cc CRPropa.h  units.h    vector3d.h \
                  Algo/algo_vec.h \
		  Exceptions/crp_err.h Exceptions/fits_err.h Exceptions/xml_err.h


