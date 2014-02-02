#ifndef __FUNCTION_H
#define __FUNCTION_H

#include <boost/python.hpp>
#include <boost/function.hpp>
#include <math.h>

namespace bp = boost::python;

class Function {

public:

  int CONST_FUNC;
  double TC, EC, Alpha;

  boost::function<double (double)> TransientVoltage;

  // Constructor
  Function();

  // Functions
  double Eval_TC(double T);
  double Eval_dTCdT(double T);
  double Eval_EC(double T);
  double Eval_Alpha(double T);
  void SetTransientVoltage(bp::object func);
  static double DefaultTransientVoltage(double t);

};

class FunctionWrapper {

  bp::object _callable;

public:

  FunctionWrapper(bp::object function) : _callable(function) {};

  double operator()(double t);

};
#endif
