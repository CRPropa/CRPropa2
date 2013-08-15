/**
	@file       defqueue.h
	@author     Tristan Beau, beau@in2p3.fr
	@brief      Container definition for QUEUE<TParticle>

*/

#ifndef _DEFQUEUE_H_
#define _DEFQUEUE_H_

#define QUEUE deque
//#define QUEUE queue // container adaptator
//#define QUEUE vector 

#if QUEUE==deque
#include <deque>

#elif QUEUE==queue
#include <queue>

#elif QUEUE==vector
#include <vector>

#endif

#endif
