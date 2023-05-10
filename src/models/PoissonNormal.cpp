#include <TMB.hpp>				// Links in the TMB libraries

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(y);				                // Data vector transmitted from R
  DATA_VECTOR(x);                       // Data vector transmitted from R
  DATA_MATRIX(X);                       // Design matrix transmitted from R
  
  PARAMETER_VECTOR(u);			            // Random effects
  
  // Parameters
  PARAMETER_VECTOR(beta);         // Parameter value transmitted from R
  PARAMETER(log_sigma_u);				// Parameter value transmitted from R
  
  vector<Type> lambda  = exp(X*beta-log(x)+u);
  Type sigma_u = exp(log_sigma_u);
  
  int nobs = y.size();
  Type mean_ran = Type(0);
  
  Type f = 0;                           // Declare the "objective function"
  for(int t=0; t < nobs; t++){
    f -= dnorm(u[t],mean_ran,sigma_u,true);
    f -= dpois(y[t],lambda[t],true);
  }
  
  return f;
}
