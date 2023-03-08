#include <TMB.hpp>				// Links in the TMB libraries

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(y);				                // Data vector transmitted from R
  DATA_VECTOR(w)                        // Data vector transmitted from R
  DATA_FACTOR(ageGroup);                // Data factor transmitted from R

  PARAMETER_VECTOR(u);			            // Random effects
   
  // Parameters
  PARAMETER_VECTOR(log_lambda);         // Parameter value transmitted from R
  PARAMETER_VECTOR(log_sigma_u);				// Parameter value transmitted from R
  
  vector<Type> lambda  = exp(log_lambda);
  vector<Type> sigma_u = exp(log_sigma_u);

  int nobs = y.size();
  Type mean_ran = Type(0);
  
  int i;

  Type f = 0;                           // Declare the "objective function"
  for(int t=0; t < nobs; t++){
    i = ageGroup[t];
    f -= dnorm(u[t],mean_ran,sigma_u[i],true);
    f -= dpois(y[t],exp(lambda[i]-log(w[t]))*exp(u[t]),true);
  }
  
  return f;
}
