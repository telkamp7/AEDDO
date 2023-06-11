#include <TMB.hpp>				// Links in the TMB libraries

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(y);				                // Data vector transmitted from R
  DATA_VECTOR(n);                       // Data vector transmitted from R
  DATA_MATRIX(X);                       // Design matrix transmitted from R
  
  // Parameters
  PARAMETER_VECTOR(beta);         // Parameter value transmitted from R
  PARAMETER(log_phi_u);				// Parameter value transmitted from R
  
  vector<Type> lambda  = exp(X*beta-log(n)); // Construct the model parameters
  Type phi_u = exp(log_phi_u); // ... and the model parameters
  
  Type r = 1/phi_u; // Construct the size
  vector<Type> p = 1/(lambda*phi_u+1); // ... and the probability parameter
  
  Type f = -sum(dnbinom(y, r, p,true)); // Calculate the "objective function"
  
  return f;
}
