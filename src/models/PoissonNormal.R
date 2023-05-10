#' 
#'
#'

PoissonNormal <- function(
    data,
    parameters,
    hessian = TRUE,
    silent = TRUE,
    method = "BFGS",
    DLL,
    cpp
){
  
  # Construct obj. function with derivatives based on a compiled C++ template.
  PoisN <- TMB::MakeADFun(
    data = data,
    parameters = parameters,
    hessian = hessian,
    DLL = libraryName,
    silent = silent,
    method = method
  )
  
  # Optimize the function
  opt <- do.call(what = "optim", args = PoisG)
  
  return(opt)
}
