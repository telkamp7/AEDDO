#include <TMB.hpp>				// Links in the TMB libraries

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(y);				    // Data vector transmitted from R
  
  PARAMETER_VECTOR(u);			// Random effects

  // Parameters
  PARAMETER(lambda);
  PARAMETER(phi);

  int nobs = y.size();			// Number of time points
  
  Type f = 0;       // Declare the "objective function" (neg. log. likelihood)
  
  for(int i=0; i < nobs; i++){
    f -= dnorm(u[i],Type(0),Type(1),true);
    Type v = pnorm(u[i],Type(0),Type(1));
    Type w = qgamma(v,Type(1),phi);
    f -= dpois(y[i],lambda*w,true);
  }
  
  return f;
}
