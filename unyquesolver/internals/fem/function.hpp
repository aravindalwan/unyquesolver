#ifndef __FUNCTION_H
#define __FUNCTION_H

#include <math.h>

class Function {

public:

  int CONST_FUNC;
  double TC, EC, Alpha;

  // Constructor
  Function();

  // Functions
  double Eval_TC(double T);
  double Eval_dTCdT(double T);
  double Eval_EC(double T);
  double Eval_Alpha(double T);

};
#endif
