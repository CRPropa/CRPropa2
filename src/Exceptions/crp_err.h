
/**
	@file		crp_err.h
	@author		Tristan Beau, beau@in2p3.fr
	@brief		Parent class for all sub class used for exceptions
*/

#ifndef _CRP_ERR_H_
#define _CRP_ERR_H_

enum MainErrors {
  ERR_MAIN = 0x1,
  ERR_CRP,
  ERR_XML,
  ERR_FITS
};

#include <iostream>
#include <string>
#include <stdlib.h>

/**
   @class TCrpErr
   @brief Basic class for exception handling
*/

class TCrpErr {
 public:
  TCrpErr() { }
  TCrpErr(const char *aMessage) { 
    std::cerr << aMessage << std::endl;
    std::cerr << "Crp exception catched. " 
	      << "It is a major one but without any management. " 
	      << "Exiting." << std::endl;
    exit(ERR_CRP);
  } /**<  Prints a message to error stream, and leave the simulation without handling */
  TCrpErr(std::string aMessage) {
    std::cerr << aMessage << std::endl;
    std::cerr << "Crp exception catched. " 
	      << "It is a major one but without any management. " 
	      << "Exiting." << std::endl;
    exit(ERR_CRP);
  } /**<  Prints a message to error stream, and leave the simulation without handling */
  ~TCrpErr() { }
} ;

// TO DO: avant le exit, fermer les fichiers ascii ou fits dans la 
//mesure du possible, histoire de sauvegarder qd meme quelque chose...


#endif

