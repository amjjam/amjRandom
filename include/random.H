/******************************************************************************
 * This defines class RANDOM which is a random number generator class         *
 ******************************************************************************/

#include <math.h>

class RANDOM{
public:
  RANDOM();
  RANDOM(long seed);
  ~RANDOM();
  float uniform();
  float gaussian();
  float normal();
private:
  long seed;
  float ran1();
  float gasdev();
};
