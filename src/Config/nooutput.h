
/**
   @file   nooutput.h
   @author Eric Armengaud, armengau@in2p3.fr
   @brief  Empty class in case no output is required
*/

#ifndef _NOOUTPUT_H_
#define _NOOUTPUT_H_

#include "outputdata.h"

class TNoOutput : public TOutputData {

 public:
	TNoOutput() { SetType(OUTPUT_NO);}
  ~TNoOutput() {}
 
};

#endif
