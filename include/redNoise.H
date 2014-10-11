/******************************************************************************
 * This file defines the class REDNOISE, which defines a red noise generator  *
 * The red noise can be a single value or an array of values, and the memory  *
 * can be provided as a pointer or can be allocated internally                *
 ******************************************************************************/

#ifndef _REDNOISE_H_
#define _REDNOISE_H_

#include <math.h>
#include "random.h"

class REDNOISE{
public:
  REDNOISE(RANDOM *random, int nVals, float tau, float sigma, float value);
  REDNOISE(RANDOM *random, int nVals, float tau, float sigma, float *values);
  ~REDNOISE();
  void advance(float dt);
  float getValue(int i);
  float *getValues();
private:
  int nVals;
  int internal;
  float *alpha;
  float *sigma;
  float *value;

  RANDOM *random;
};

#endif
