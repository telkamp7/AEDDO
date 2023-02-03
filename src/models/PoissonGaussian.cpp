#include <TMB.hpp>				// Links in the TMB libraries

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(y);				// Data vector transmitted from R
  PARAMETER_VECTOR(u);				// Random effects
  
  // Parameters
  PARAMETER(lambda);				// Parameter value transmitted from R
  PARAMETER(sigma_u);				// Parameter value transmitted from R
  
  int n = y.size();
  Type mean_ran = Type(0);
  
  Type f = 0;
  
  for(int j=0; j < n; j++){
    f -= dnorm(u[j],mean_ran,sigma_u,true);
  }
  
  for(int i=0; i < n; i++){
    f -= dpois(y[i],lambda*exp(u[i]),true);
  }
  
  return f;
}
