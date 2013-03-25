#include "function.hpp"

//------------------------------------------------------------------------------
Function::Function() {
  // Default is to allow varying functions
  CONST_FUNC = 0;

  // Initialize constant values of functions
  TC = 146.4; EC = 2381; Alpha = 2.568e-6;
}
//------------------------------------------------------------------------------
double Function::Eval_TC(double T) {

  //Initialize coefficients needed for this function
  double a = 615.3, b = -0.006284, c = 64.19, d = -0.0006563;

  // Return value of thermal conductivity at the given temperature
  // TC = a*exp(b*T)+c*exp(d*T)
  if (CONST_FUNC == 1)
    return TC;
  else
    return (a*exp(b*T)+c*exp(d*T));

}
//------------------------------------------------------------------------------
double Function::Eval_dTCdT(double T) {

  //Initialize coefficients needed for this function
  double a = 615.3, b = -0.006284, c = 64.19, d = -0.0006563;

  // Return value of gradient of thermal conductivity at the given temperature
  // dTCdT = a*b*exp(b*T)+c*d*exp(d*T)
  if (CONST_FUNC == 1)
    return 0.0;
  else
    return (a*b*exp(b*T)+c*d*exp(d*T));

}
//------------------------------------------------------------------------------
double Function::Eval_EC(double T) {

  // Return value of electrical conductivity at the given temperature
  return EC;

}
//------------------------------------------------------------------------------
double Function::Eval_Alpha(double T) {

  //Initialize coefficients needed for this function
  double a = 3.776e-6, b = 0.0001257, c = -7.627e-6, d = -0.005766;

  // Return value of electrical conductivity at the given temperature
  if (CONST_FUNC == 1)
    return Alpha;
  else
    return (a*exp(b*T)+c*exp(d*T));

}
