/*
 * Automatically Generated from Mathematica.
 * Tue 24 Oct 2023 17:04:57 GMT-04:00
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
  double t2443;
  double t2481;
  double t2368;
  double t2408;
  double t2360;
  double t2400;
  double t2445;
  double t2452;
  double t2467;
  double t2468;
  double t2469;
  double t2491;
  double t2492;
  double t2493;
  double t2497;
  double t2498;
  double t2499;
  double t2500;
  double t2501;
  double t2502;
  double t2503;
  double t2484;
  double t2485;
  double t2486;
  double t2487;
  double t2488;
  double t2489;
  double t2505;
  double t2510;
  double t2511;
  double t2512;
  double t2513;
  double t2504;
  double t2506;
  double t2514;
  double t2369;
  double t2426;
  double t2427;
  double t2428;
  double t2494;
  double t2526;
  double t2527;
  double t2528;
  double t2550;
  double t2551;
  double t2552;
  double t2542;
  double t2543;
  double t2544;
  double t2546;
  double t2547;
  double t2548;
  double t2567;
  double t2568;
  double t2569;
  double t2571;
  double t2572;
  double t2573;
  double t2490;
  double t2495;
  double t2525;
  double t2529;
  double t2530;
  double t2531;
  double t2534;
  double t2535;
  double t2536;
  double t2537;
  double t2538;
  double t2539;
  double t2549;
  double t2553;
  double t2597;
  double t2598;
  double t2558;
  double t2600;
  double t2601;
  double t2560;
  t2443 = Cos(var1[6]);
  t2481 = Sin(var1[6]);
  t2368 = Sin(var1[2]);
  t2408 = Sin(var1[5]);
  t2360 = Cos(var1[5]);
  t2400 = Cos(var1[2]);
  t2445 = -1.*t2443;
  t2452 = 1. + t2445;
  t2467 = 0.4*t2452;
  t2468 = 0.64*t2443;
  t2469 = t2467 + t2468;
  t2491 = t2360*t2443;
  t2492 = -1.*t2408*t2481;
  t2493 = t2491 + t2492;
  t2497 = t2469*t2443;
  t2498 = Power(t2481,2);
  t2499 = 0.24*t2498;
  t2500 = t2497 + t2499;
  t2501 = -1.*t2443*t2408;
  t2502 = -1.*t2360*t2481;
  t2503 = t2501 + t2502;
  t2484 = t2469*t2481;
  t2485 = -0.24*t2443*t2481;
  t2486 = t2484 + t2485;
  t2487 = t2443*t2408;
  t2488 = t2360*t2481;
  t2489 = t2487 + t2488;
  t2505 = -1.*t2368*t2493;
  t2510 = -1.*t2400*t2360;
  t2511 = t2368*t2408;
  t2512 = t2510 + t2511;
  t2513 = -0.11*t2512;
  t2504 = t2400*t2503;
  t2506 = t2504 + t2505;
  t2514 = -1.*t2368*t2503;
  t2369 = -1.*t2360*t2368;
  t2426 = -1.*t2400*t2408;
  t2427 = t2369 + t2426;
  t2428 = -0.11*t2427;
  t2494 = t2400*t2493;
  t2526 = -1.*t2360*t2443;
  t2527 = t2408*t2481;
  t2528 = t2526 + t2527;
  t2550 = t2360*t2469;
  t2551 = -0.24*t2408*t2481;
  t2552 = t2550 + t2551;
  t2542 = -1.*t2469*t2408;
  t2543 = -0.24*t2360*t2481;
  t2544 = t2542 + t2543;
  t2546 = t2469*t2408;
  t2547 = 0.24*t2360*t2481;
  t2548 = t2546 + t2547;
  t2567 = -1.*t2469*t2481;
  t2568 = 0.24*t2443*t2481;
  t2569 = t2567 + t2568;
  t2571 = Power(t2443,2);
  t2572 = -0.24*t2571;
  t2573 = t2497 + t2572;
  t2490 = -1.*t2368*t2489;
  t2495 = t2490 + t2494;
  t2525 = -1.*t2486*t2506;
  t2529 = t2400*t2528;
  t2530 = t2514 + t2529;
  t2531 = -1.*t2500*t2530;
  t2534 = t2368*t2503;
  t2535 = t2534 + t2494;
  t2536 = -1.*t2486*t2535;
  t2537 = t2368*t2528;
  t2538 = t2504 + t2537;
  t2539 = -1.*t2500*t2538;
  t2549 = t2548*t2493;
  t2553 = t2503*t2552;
  t2597 = -0.24*t2443*t2408;
  t2598 = t2597 + t2543;
  t2558 = -1.*t2503*t2548;
  t2600 = 0.24*t2360*t2443;
  t2601 = t2600 + t2551;
  t2560 = -1.*t2552*t2528;
  p_output1[0]=0;
  p_output1[1]=0;
  p_output1[2]=(-0.5*(t2428 - 1.*t2486*t2495 - 1.*t2500*t2506)*var2[0] - 0.5*(-1.*t2486*(-1.*t2400*t2489 + t2505) + t2513 - 1.*t2500*(-1.*t2400*t2493 + t2514))*var2[1])*var2[5];
  p_output1[3]=0;
  p_output1[4]=0;
  p_output1[5]=(-0.5*(t2428 + t2536 + t2539)*var2[0] - 0.5*(t2513 + t2525 + t2531)*var2[1] - 0.5*(-1.*t2500*(t2493*t2544 + t2549 + t2489*t2552 + t2553) - 1.*t2486*(-1.*t2503*t2544 - 1.*t2493*t2552 + t2558 + t2560))*var2[2])*var2[5];
  p_output1[6]=var2[5]*(-0.5*(t2536 + t2539 - 1.*t2535*t2569 - 1.*(t2400*t2489 + t2368*t2493)*t2573)*var2[0] - 0.5*(t2525 + t2531 - 1.*t2506*t2569 - 1.*t2495*t2573)*var2[1] - 0.5*(-1.*(t2489*t2548 + t2493*t2552)*t2569 - 1.*(-1.*t2493*t2548 - 1.*t2503*t2552)*t2573 - 1.*t2500*(t2549 + t2553 + t2493*t2598 + t2489*t2601) - 1.*t2486*(t2558 + t2560 - 1.*t2503*t2598 - 1.*t2493*t2601))*var2[2] - 0.5*(-2.*t2500*t2569 - 2.*t2486*t2573)*var2[5] + 0.12*t2569*var2[6]);
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

#include "Ce3_vec6_five_link_walker.hh"

namespace SymFunction
{

void Ce3_vec6_five_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
