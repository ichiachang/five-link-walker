/*
 * Automatically Generated from Mathematica.
 * Tue 24 Oct 2023 17:04:55 GMT-04:00
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
  double t2348;
  double t2333;
  double t2334;
  double t2349;
  double t2353;
  double t2335;
  double t2350;
  double t2351;
  double t2309;
  double t2354;
  double t2358;
  double t2359;
  double t2429;
  double t2438;
  double t2440;
  double t2441;
  double t2442;
  double t2460;
  double t2461;
  double t2462;
  double t2352;
  double t2388;
  double t2470;
  double t2471;
  double t2472;
  double t2444;
  double t2456;
  double t2457;
  double t2458;
  double t2459;
  double t2463;
  double t2464;
  double t2465;
  double t2466;
  double t2473;
  double t2474;
  double t2475;
  double t2476;
  double t2477;
  double t2478;
  t2348 = Cos(var1[3]);
  t2333 = Cos(var1[4]);
  t2334 = Sin(var1[3]);
  t2349 = Sin(var1[4]);
  t2353 = Cos(var1[2]);
  t2335 = -1.*t2333*t2334;
  t2350 = -1.*t2348*t2349;
  t2351 = t2335 + t2350;
  t2309 = Sin(var1[2]);
  t2354 = t2348*t2333;
  t2358 = -1.*t2334*t2349;
  t2359 = t2354 + t2358;
  t2429 = -1.*t2333;
  t2438 = 1. + t2429;
  t2440 = 0.4*t2438;
  t2441 = 0.64*t2333;
  t2442 = t2440 + t2441;
  t2460 = t2348*t2442;
  t2461 = -0.24*t2334*t2349;
  t2462 = t2460 + t2461;
  t2352 = -1.*t2309*t2351;
  t2388 = t2353*t2351;
  t2470 = -1.*t2348*t2333;
  t2471 = t2334*t2349;
  t2472 = t2470 + t2471;
  t2444 = -0.24*t2348*t2349;
  t2456 = t2442*t2334;
  t2457 = 0.24*t2348*t2349;
  t2458 = t2456 + t2457;
  t2459 = t2458*t2359;
  t2463 = t2351*t2462;
  t2464 = t2333*t2334;
  t2465 = t2348*t2349;
  t2466 = t2464 + t2465;
  t2473 = t2353*t2472;
  t2474 = t2352 + t2473;
  t2475 = 0.12*var2[1]*t2474;
  t2476 = t2309*t2472;
  t2477 = t2388 + t2476;
  t2478 = 0.12*var2[0]*t2477;
  p_output1[0]=0;
  p_output1[1]=0;
  p_output1[2]=(0.12*(-1.*t2309*t2359 + t2388)*var2[0] + 0.12*(t2352 - 1.*t2353*t2359)*var2[1])*var2[4];
  p_output1[3]=(t2475 + t2478 + 0.12*(t2359*(-1.*t2334*t2442 + t2444) + t2459 + t2463 + t2462*t2466)*var2[2])*var2[4];
  p_output1[4]=(t2475 + t2478 + 0.12*(t2359*(-0.24*t2333*t2334 + t2444) + t2459 + t2463 + (0.24*t2333*t2348 + t2461)*t2466)*var2[2] + 0.12*(0.24*t2333*t2349 - 1.*t2349*t2442)*var2[3])*var2[4];
  p_output1[5]=0;
  p_output1[6]=0;
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

#include "Ce3_vec5_five_link_walker.hh"

namespace SymFunction
{

void Ce3_vec5_five_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
