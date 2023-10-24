/*
 * Automatically Generated from Mathematica.
 * Tue 24 Oct 2023 17:04:17 GMT-04:00
 */

#ifdef MATLAB_MEX_FILE
#include <stdexcept>
#include <cmath>
#include<math.h>
/**
 * Copied from Wolfram Mathematica C Definitions file mdefs.hpp
 * Changed marcos to inline functions (Eric Cousineau)
 */
inline double Power(double x, double y) { return pow(x, y); }
inline double Sqrt(double x) { return sqrt(x); }

inline double Abs(double x) { return fabs(x); }

inline double Exp(double x) { return exp(x); }
inline double Log(double x) { return log(x); }

inline double Sin(double x) { return sin(x); }
inline double Cos(double x) { return cos(x); }
inline double Tan(double x) { return tan(x); }

inline double ArcSin(double x) { return asin(x); }
inline double ArcCos(double x) { return acos(x); }
inline double ArcTan(double x) { return atan(x); }

/* update ArcTan function to use atan2 instead. */
inline double ArcTan(double x, double y) { return atan2(y,x); }

inline double Sinh(double x) { return sinh(x); }
inline double Cosh(double x) { return cosh(x); }
inline double Tanh(double x) { return tanh(x); }

const double E	= 2.71828182845904523536029;
const double Pi = 3.14159265358979323846264;
const double Degree = 0.01745329251994329576924;

inline double Sec(double x) { return 1/cos(x); }
inline double Csc(double x) { return 1/sin(x); }

#endif

/*
 * Sub functions
 */
static void output1(double *p_output1,const double *var1)
{
  double t4;
  double t8;
  double t9;
  double t10;
  double t11;
  double t12;
  double t13;
  double t15;
  double t19;
  double t20;
  double t21;
  double t22;
  double t28;
  double t29;
  double t30;
  double t31;
  double t32;
  double t34;
  double t38;
  double t39;
  double t40;
  double t41;
  double t16;
  double t17;
  double t55;
  double t56;
  double t57;
  double t35;
  double t36;
  double t65;
  double t66;
  double t67;
  t4 = Sin(var1[2]);
  t8 = Cos(var1[3]);
  t9 = t8*t4;
  t10 = Cos(var1[2]);
  t11 = Sin(var1[3]);
  t12 = t10*t11;
  t13 = t9 + t12;
  t15 = Cos(var1[4]);
  t19 = t10*t8;
  t20 = -1.*t4*t11;
  t21 = t19 + t20;
  t22 = Sin(var1[4]);
  t28 = Cos(var1[5]);
  t29 = t28*t4;
  t30 = Sin(var1[5]);
  t31 = t10*t30;
  t32 = t29 + t31;
  t34 = Cos(var1[6]);
  t38 = t10*t28;
  t39 = -1.*t4*t30;
  t40 = t38 + t39;
  t41 = Sin(var1[6]);
  t16 = -1.*t15;
  t17 = 1. + t16;
  t55 = -1.*t8*t4;
  t56 = -1.*t10*t11;
  t57 = t55 + t56;
  t35 = -1.*t34;
  t36 = 1. + t35;
  t65 = -1.*t28*t4;
  t66 = -1.*t10*t30;
  t67 = t65 + t66;
  p_output1[0]=0.03125*(0.11*t13 + 0.4*t13*t17 - 0.4*t21*t22 + 0.64*(t13*t15 + t21*t22) + 0.11*t32 + 0.4*t32*t36 - 0.4*t40*t41 + 0.64*(t32*t34 + t40*t41) + 4.*var1[0] + 28.*(0.24*t4 + var1[0]));
  p_output1[1]=0;
  p_output1[2]=0.03125*(0.11*t21 + 0.4*t17*t21 + 0.11*t40 + 0.4*t36*t40 - 0.4*t22*t57 + 0.64*(t15*t21 + t22*t57) - 0.4*t41*t67 + 0.64*(t34*t40 + t41*t67) + 4.*var1[1] + 28.*(0.24*t10 + var1[1]));
}



#ifdef MATLAB_MEX_FILE

#include "mex.h"
/*
 * Main function
 */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
  size_t mrows, ncols;

  double *var1;
  double *p_output1;

  /*  Check for proper number of arguments.  */ 
  if( nrhs != 1)
    {
      mexErrMsgIdAndTxt("MATLAB:MShaped:invalidNumInputs", "One input(s) required (var1).");
    }
  else if( nlhs > 1)
    {
      mexErrMsgIdAndTxt("MATLAB:MShaped:maxlhs", "Too many output arguments.");
    }

  /*  The input must be a noncomplex double vector or scaler.  */
  mrows = mxGetM(prhs[0]);
  ncols = mxGetN(prhs[0]);
  if( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) ||
    ( !(mrows == 7 && ncols == 1) && 
      !(mrows == 1 && ncols == 7))) 
    {
      mexErrMsgIdAndTxt( "MATLAB:MShaped:inputNotRealVector", "var1 is wrong.");
    }

  /*  Assign pointers to each input.  */
  var1 = mxGetPr(prhs[0]);
   


   
  /*  Create matrices for return arguments.  */
  plhs[0] = mxCreateDoubleMatrix((mwSize) 1, (mwSize) 3, mxREAL);
  p_output1 = mxGetPr(plhs[0]);


  /* Call the calculation subroutine. */
  output1(p_output1,var1);


}

#else // MATLAB_MEX_FILE

#include "pos_CoM_func.hh"

namespace SymFunction
{

void pos_CoM_func_raw(double *p_output1, const double *var1)
{
  // Call Subroutines
  output1(p_output1, var1);

}

}

#endif // MATLAB_MEX_FILE
