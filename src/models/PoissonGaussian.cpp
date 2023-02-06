#include <TMB.hpp>				// Links in the TMB libraries

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(y);				        // Data vector transmitted from R
  DATA_FACTOR(ageLabel);        // Data factor transmitted from R

  PARAMETER_VECTOR(u);			    // Random effects
   
  // Parameters
  PARAMETER_VECTOR(lambda); 	  // Parameter value transmitted from R
  PARAMETER(log_sigma_u);				// Parameter value transmitted from R
  
  Type sigma_u = exp(log_sigma_u);

  ADREPORT(sigma_u);

  int nobs = y.size();
  Type mean_ran = Type(0);
  
  int j;

  Type f = 0;               // Declare the "objective function" (neg. log. likelihood)
  for(int i=0; i < nobs; i++){
    f -= dnorm(u[i],mean_ran,sigma_u,true);
    j = ageLabel[i];
    f -= dpois(y[i],lambda[j]*exp(u[i]),true);
  }
  
  return f;
}
