#include <TMB.hpp>				// Links in the TMB libraries

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(y);				// Data vector transmitted from R
  PARAMETER_VECTOR(u);				// Random effects
  
  // Parameters
  PARAMETER(lambda);
  //PARAMETER(alpha);
  PARAMETER(phi);
  //PARAMETER_VECTOR(ageLabel);
  //PARAMETER_VECTOR(landsdel);
  
  Type alpha = Type(1);
  int n = y.size();				// Number of time points
  
  Type f =  0;					// Declare the "objective function" (neg. log. likelihood)
  for(int t = 1; t < n; t++){			// start at t = 1
    f -= dgamma(u[t],alpha,phi, true);
    f -= dpois(y[t],lambda*u[t],true);
  }
  
  return f;
}
