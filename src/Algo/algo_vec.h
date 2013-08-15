
/**
	@file       algo_vec.h
	@author     Tristan Beau, beau@in2p3.fr
	@brief      Some useful macro definitions for stl-vector
                    manipulations

*/

#ifndef _ALGO_VEC_H_
#define _ALGO_VEC_H_

#define A_EQ_B_PLUS_C(AA,BB,CC) \
  transform ( (BB).begin(), \
              (BB).end(),   \
              (CC).begin(), \
              (AA).begin(), \
              plus<double>() ) 

#define A_EQ_B_MULT_K(AA,BB,KK) \
  transform ( (BB).begin(), \
              (BB).end(),   \
              (AA).begin(), \
              bind1st(multiplies<double>(), (KK) ) )

#endif

