/******************************************************************************
 * This defines class RANDOM which is a random number generator class         *
 ******************************************************************************/

#include "../include/random.H"

/*=============================================================================
  RANDOM::RANDOM() - creator. Creates the random number generator with
  the default initial seed of -1.
  ============================================================================*/
RANDOM::RANDOM(){
  seed=1;
  iy=0;
  iset=0;
}


/*=============================================================================
  RANDOM::RANDOM(long seed) - creator. Creates the random number
  generator with an initial seed.

  long seed - the initial seed. Must be a negative number.
  ============================================================================*/
RANDOM::RANDOM(long seed){
  RANDOM::seed=seed;
  iy=0;
  iset=0;
}


/*=============================================================================
  float RANDOM::uniform() - returns a random number uniformly
  distributed in the interval [0;1[
  ============================================================================*/
float RANDOM::uniform(){
  return ran1();
}


/*=============================================================================
  float RANDOM:gaussian() - returns a random number normally/gaussian
  distributed with a mean of zero and a standard deviation of one. To
  create a random number with a different offset and standard
  deviation, multiply by the standard deviation and then add the
  offset.
  ============================================================================*/
float RANDOM::gaussian(){
  return gasdev();
}


/*=============================================================================
  float RANDOM::normal() - synonym for gaussian()
  ============================================================================*/
float RANDOM::normal(){
  return gaussian();
}


/******************************************************************************
 ******************************************************************************
 *                      Private members                                       *
 ******************************************************************************
 ******************************************************************************/


/*=============================================================================
  float RANDOM::ran1() - return a number uniformly distributed in the
  [0;1[ interval. This function is from Numerical Recipes.
  ============================================================================*/
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)
float RANDOM::ran1(){
  int j;
  long k;
  /*  static long iy=0;
      static long iv[NTAB];*/
  float temp;
  
  if (seed <= 0 || !iy) {
    if (-seed < 1) seed=1;
    else seed = -seed;
    for (j=NTAB+7;j>=0;j--) {
      k=seed/IQ;
      seed=IA*(seed-k*IQ)-IR*k;
      if (seed < 0) seed += IM;
      if (j < NTAB) iv[j] = seed;
    }
    iy=iv[0];
  }
  k=seed/IQ;
  seed=IA*(seed-k*IQ)-IR*k;
  if (seed < 0) seed += IM;
  j=iy/NDIV;
  iy=iv[j];
  iv[j] = seed;
  if ((temp=AM*iy) > RNMX) return RNMX;
  else return temp;
}
#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX


/*=============================================================================
  float RANDOM::gasdev() - returns a random number from a normal
  distribution with zero mean and unity standard deviation.
  ============================================================================*/
float RANDOM::gasdev(){
  /*static int iset=0;
    static float gset;*/
  float fac,rsq,v1,v2;
  
  if  (iset == 0) {
    do {
      v1=2.0*ran1()-1.0;
      v2=2.0*ran1()-1.0;
      rsq=v1*v1+v2*v2;
    } while (rsq >= 1.0 || rsq == 0.0);
    fac=sqrt(-2.0*log(rsq)/rsq);
    gset=v1*fac;
    iset=1;
    return v2*fac;
  } else {
    iset=0;
    return gset;
  }
}
