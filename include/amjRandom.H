/******************************************************************************
 * This defines class amjRandom which is a random number generator class      *
 ******************************************************************************/

#ifndef _AMJRANDOM_H_
#define _AMJRANDOM_H_

#include <math.h>

class amjRandom{
public:
  amjRandom();
  amjRandom(long seed);
  ~amjRandom();
  float uniform();
  float gaussian();
  float normal();
  float poisson(float x);
private:
  long seed;
  long iy;
  long iv[32];
  int iset;
  float gset;
  float ran1();
  float gammln(float xx);
  float gasdev();
  float poidev(float xm);
};

#endif
