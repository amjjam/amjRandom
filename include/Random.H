/******************************************************************************
 * This defines class RANDOM which is a random number generator class         *
 ******************************************************************************/

#ifndef _RANDOM_H_
#define _RANDOM_H_

#include <math.h>

class Random{
public:
  Random();
  Random(long seed);
  ~Random();
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
