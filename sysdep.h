
#ifndef SYS_DEP_H
#define SYS_DEP_H

#if HAVE_CONFIG_H

  #ifdef PACKAGE
    #undef PACKAGE
  #endif

  #ifdef PACKAGE_BUGREPORT
    #undef PACKAGE_BUGREPORT
  #endif 

  #ifdef PACKAGE_VERSION
    #undef PACKAGE_VERSION
  #endif

  #ifdef PACKAGE_NAME
	  #undef PACKAGE_NAME
	#endif

  #ifdef PACKAGE_STRING
    #undef PACKAGE_STRING
	#endif

  #ifdef PACKAGE_TARNAME
	  #undef PACKAGE_TARNAME 
	#endif

  #ifdef VERSION
    #undef VERSION
  #endif

  #include "config.h"

#endif //HAVE_CONFIG_H

#endif

