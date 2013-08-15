/**
   @file    nullmagfield.h
   @author  Eric Armengaud, armengau@in2p3.fr
   @brief   Trivial class describing a null magnetic field
*/

#ifndef _NULLMAGFIELD_H_
#define _NULLMAGFIELD_H_

#include "magfield.h"

/**
   @class TNullMagField
   @brief No magnetic field at all.
*/


class TNullMagField : public TMagField {
 public:
  TNullMagField() { SetType(MAGFIELD_NO); }
  ~TNullMagField() {}
};


#endif
