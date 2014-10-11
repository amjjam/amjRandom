/******************************************************************************
 * This file defines the class REDNOISE, which defines a red noise generator  *
 * The red noise can be a single value or an array of values, and the memory  *
 * can be provided as a pointer or can be allocated internally                *
 ******************************************************************************/

#include "../include/redNoise.H"

/*=============================================================================
  REDNOISE::REDNOISE(RANDOM *random, int nVals, float tau, float
  sigma, float value) - constructor for red noise generator with nVals
  time-series stored internally.

  RANDOM *random - pointer to the random number generator
  int nVals - number of red noise time-series to generate
  float tau - the time-constant over which one sigma change occurs
  float sigma - the change over one time-constant
  float value - the initial value of the random time-series
  ============================================================================*/
REDNOISE::REDNOISE(RANDOM *random, int nVals, float tau, float sigma, 
		   float value){
  
  REDNOISE::random=random;
  REDNOISE::nVals=nVals;
  internal=1;
  REDNOISE::alpha=new float[nVals];
  REDNOISE::sigma=new float[nVals];
  REDNOISE::value=new float[nVals];

  for(int i=0;i<nVals;i++){
    REDNOISE::alpha[i]=1.0-1.0/tau;
    REDNOISE::sigma[i]=sigma;
    REDNOISE::value[i]=value;
  }
}

/*=============================================================================
  REDNOISE::REDNOISE(RANDOM *random, int nVals, float tau, float
  sigma, float *values) - constructor for red noise generator with
  nVals time-series with the values provided as a pointer

  RANDOM *random - pointer to the random number generator
  int nVals - number of red noise time-series to generate
  float tau - the time-constant over which one sigma change occurs
  float sigma - the change over one time-constant
  float *valueS - the initial value of the random time-series
  ============================================================================*/
REDNOISE::REDNOISE(RANDOM *random, int nVals, float tau, float sigma, 
		   float *values){
  
  REDNOISE::random=random;
  REDNOISE::nVals=nVals;
  internal=0;
  REDNOISE::alpha=new float[nVals];
  REDNOISE::sigma=new float[nVals];
  REDNOISE::value=new float[nVals];
  
  for(int i=0;i<nVals;i++){
    REDNOISE::alpha[i]=1.0-1.0/tau;
    REDNOISE::sigma[i]=sigma;
  }
  REDNOISE::value=values;
}


/*=============================================================================
  REDNOISE::~REDNOISE() - destructor. Deacllocates memory
  ============================================================================*/
REDNOISE::~REDNOISE(){
  delete [] alpha;
  delete [] sigma;
  if(internal==1) delete [] value;
}


/*=============================================================================
  float REDNOISE::advance(float dt) - advance the nVals time-series by time dt
  
  float dt - time to update forward
  ============================================================================*/
void REDNOISE::advance(float dt){
  float n=1/dt,rho;
  for(int i=0;i<nVals;i++){
    rho=sqrt((1-alpha[i])*(1-alpha[i])/
	     (n-2*alpha[i]-n*alpha[i]*alpha[i]+2*powf(alpha[i],n+1))/dt);
    value[i]+=sqrt(dt)*sigma[i]*rho*random->gaussian();
  }
}


/*=============================================================================
  float REDNOISE::getValue(int i) - returns the value of one element
  of the noise array.

  int i - which array element to return
  ============================================================================*/
float REDNOISE::getValue(int i){
  return value[i];
}


/*=============================================================================
  float *getValues() - returns pointer to the noise values array
  
  ============================================================================*/
float *REDNOISE::getValues(){

  return value;
}
