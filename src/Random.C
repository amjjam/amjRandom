/******************************************************************************
 * This defines class Random which is a random number generator class         *
 ******************************************************************************/

#include "../include/Random.H"

/*=============================================================================
  Random::Random() - creator. Creates the random number generator with
  the default initial seed of -1.
  ============================================================================*/
Random::Random(){
  seed=1;
  iy=0;
  iset=0;
}


/*=============================================================================
  Random::Random(long seed) - creator. Creates the random number
  generator with an initial seed.

  long seed - the initial seed. Must be a negative number.
  ============================================================================*/
Random::Random(long seed){
  Random::seed=seed;
  iy=0;
  iset=0;
}


/*=============================================================================
  ~Random() - destructor
  ============================================================================*/
Random::~Random(){

}


/*=============================================================================
  float Random::uniform() - returns a random number uniformly
  distributed in the interval [0;1[
  ============================================================================*/
float Random::uniform(){
  return ran1();
}


/*=============================================================================
  float Random:gaussian() - returns a random number normally/gaussian
  distributed with a mean of zero and a standard deviation of one. To
  create a random number with a different offset and standard
  deviation, multiply by the standard deviation and then add the
  offset.
  ============================================================================*/
float Random::gaussian(){
  return gasdev();
}


/*=============================================================================
  float Random::normal() - synonym for gaussian()
  ============================================================================*/
float Random::normal(){
  return gaussian();
}


/*=============================================================================
  float Random::poisson(float x) - returns a poisson random deviate
  with expectation value x
  ============================================================================*/
float Random::poisson(float x){
  return poidev(x);
}


/******************************************************************************
 ******************************************************************************
 *                      Private members                                       *
 ******************************************************************************
 ******************************************************************************/


/*=============================================================================
  float Random::ran1() - return a number uniformly distributed in the
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
float Random::ran1(){
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
  Random::gammln(float xx) - computes a gamma function. Needed by poidev
  ============================================================================*/
float Random::gammln(float xx)
{
	double x,y,tmp,ser;
	static double cof[6]={76.18009172947146,-86.50532032941677,
		24.01409824083091,-1.231739572450155,
		0.1208650973866179e-2,-0.5395239384953e-5};
	int j;

	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<=5;j++) ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}


/*=============================================================================
  float Random::gasdev() - returns a random number from a normal
  distribution with zero mean and unity standard deviation.
  ============================================================================*/
float Random::gasdev(){
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


/*=============================================================================
  float Random::poidev(float xm) - returns a poisson random deviate
  with a expectation value of xm
  ============================================================================*/
#define PI 3.141592654
float Random::poidev(float xm){
	static float sq,alxm,g,oldm=(-1.0);
	float em,t,y;

	if (xm < 12.0) {
		if (xm != oldm) {
			oldm=xm;
			g=exp(-xm);
		}
		em = -1;
		t=1.0;
		do {
			++em;
			t *= ran1();
		} while (t > g);
	} else {
		if (xm != oldm) {
			oldm=xm;
			sq=sqrt(2.0*xm);
			alxm=log(xm);
			g=xm*alxm-gammln(xm+1.0);
		}
		do {
			do {
				y=tan(PI*ran1());
				em=sq*y+xm;
			} while (em < 0.0);
			em=floor(em);
			t=0.9*(1.0+y*y)*exp(em*alxm-gammln(em+1.0)-g);
		} while (ran1() > t);
	}
	return em;
}
#undef PI
/* (C) Copr. 1986-92 Numerical Recipes Software 1)0. */


