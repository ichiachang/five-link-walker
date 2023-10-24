/*
 * Automatically Generated from Mathematica.
 * Tue 24 Oct 2023 17:04:58 GMT-04:00
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
static void output1(double *p_output1,const double *var1,const double *var2)
{
  double t2515;
  double t2507;
  double t2508;
  double t2516;
  double t2520;
  double t2509;
  double t2517;
  double t2518;
  double t2496;
  double t2521;
  double t2522;
  double t2523;
  double t2557;
  double t2559;
  double t2561;
  double t2562;
  double t2563;
  double t2578;
  double t2579;
  double t2580;
  double t2519;
  double t2540;
  double t2588;
  double t2589;
  double t2590;
  double t2565;
  double t2574;
  double t2575;
  double t2576;
  double t2577;
  double t2581;
  double t2582;
  double t2583;
  double t2584;
  double t2591;
  double t2592;
  double t2593;
  double t2594;
  double t2595;
  double t2596;
  t2515 = Cos(var1[5]);
  t2507 = Cos(var1[6]);
  t2508 = Sin(var1[5]);
  t2516 = Sin(var1[6]);
  t2520 = Cos(var1[2]);
  t2509 = -1.*t2507*t2508;
  t2517 = -1.*t2515*t2516;
  t2518 = t2509 + t2517;
  t2496 = Sin(var1[2]);
  t2521 = t2515*t2507;
  t2522 = -1.*t2508*t2516;
  t2523 = t2521 + t2522;
  t2557 = -1.*t2507;
  t2559 = 1. + t2557;
  t2561 = 0.4*t2559;
  t2562 = 0.64*t2507;
  t2563 = t2561 + t2562;
  t2578 = t2515*t2563;
  t2579 = -0.24*t2508*t2516;
  t2580 = t2578 + t2579;
  t2519 = -1.*t2496*t2518;
  t2540 = t2520*t2518;
  t2588 = -1.*t2515*t2507;
  t2589 = t2508*t2516;
  t2590 = t2588 + t2589;
  t2565 = -0.24*t2515*t2516;
  t2574 = t2563*t2508;
  t2575 = 0.24*t2515*t2516;
  t2576 = t2574 + t2575;
  t2577 = t2576*t2523;
  t2581 = t2518*t2580;
  t2582 = t2507*t2508;
  t2583 = t2515*t2516;
  t2584 = t2582 + t2583;
  t2591 = t2520*t2590;
  t2592 = t2519 + t2591;
  t2593 = 0.12*var2[1]*t2592;
  t2594 = t2496*t2590;
  t2595 = t2540 + t2594;
  t2596 = 0.12*var2[0]*t2595;
  p_output1[0]=0;
  p_output1[1]=0;
  p_output1[2]=(0.12*(-1.*t2496*t2523 + t2540)*var2[0] + 0.12*(t2519 - 1.*t2520*t2523)*var2[1])*var2[6];
  p_output1[3]=0;
  p_output1[4]=0;
  p_output1[5]=(t2593 + t2596 + 0.12*(t2523*(-1.*t2508*t2563 + t2565) + t2577 + t2581 + t2580*t2584)*var2[2])*var2[6];
  p_output1[6]=(t2593 + t2596 + 0.12*(t2523*(-0.24*t2507*t2508 + t2565) + t2577 + t2581 + (0.24*t2507*t2515 + t2579)*t2584)*var2[2] + 0.12*(0.24*t2507*t2516 - 1.*t2516*t2563)*var2[5])*var2[6];
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

  double *var1,*var2;
  double *p_output1;

  /*  Check for proper number of arguments.  */ 
  if( nrhs != 2)
    {
      mexErrMsgIdAndTxt("MATLAB:MShaped:invalidNumInputs", "Two input(s) required (var1,var2).");
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
  mrows = mxGetM(prhs[1]);
  ncols = mxGetN(prhs[1]);
  if( !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) ||
    ( !(mrows == 7 && ncols == 1) && 
      !(mrows == 1 && ncols == 7))) 
    {
      mexErrMsgIdAndTxt( "MATLAB:MShaped:inputNotRealVector", "var2 is wrong.");
    }

  /*  Assign pointers to each input.  */
  var1 = mxGetPr(prhs[0]);
  var2 = mxGetPr(prhs[1]);
   


   
  /*  Create matrices for return arguments.  */
  plhs[0] = mxCreateDoubleMatrix((mwSize) 7, (mwSize) 1, mxREAL);
  p_output1 = mxGetPr(plhs[0]);


  /* Call the calculation subroutine. */
  output1(p_output1,var1,var2);


}

#else // MATLAB_MEX_FILE

#include "Ce3_vec7_five_link_walker.hh"

namespace SymFunction
{

void Ce3_vec7_five_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
