#include <TMB.hpp>				// Links in the TMB libraries

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(y);				    // Data vector transmitted from R
  PARAMETER_VECTOR(u);			// Random effects
  
  // Parameters
  PARAMETER(lambda);
  PARAMETER(phi);

  int nobs = y.size();				  // Number of time points
  
  Type f =  Type(0.0);					    // Declare the "objective function" (neg. log. likelihood)
  f -= sum(dnorm(u,Type(0),Type(1),true));
  //f -= dnorm(u,Type(0),Type(1),true).sum();       // Assign N(0,1) distribution u 
  vector<Type> v = pnorm(u,Type(0),Type(1));  // Uniformly distributed variables (on [0,1])
  vector<Type> w = qgamma(v,Type(1),phi);
  f -= dpois(y,lambda*w,true).sum();
  //for(int t = 1; t < nobs; t++){			// start at t = 1
    //f -= -(log(pow(phi,alpha))+log(pow(u[t],(alpha-1)))-phi*u[t]-lgamma(vector<Type>(alpha)));
    //f -= dpois(y[t],lambda*w[t],true);
  //}
  
  return f;
}
