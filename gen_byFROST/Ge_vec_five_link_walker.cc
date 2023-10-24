/*
 * Automatically Generated from Mathematica.
 * Tue 24 Oct 2023 17:05:00 GMT-04:00
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
  double t2524;
  double t2533;
  double t2541;
  double t2545;
  double t2554;
  double t2555;
  double t2556;
  double t2566;
  double t2587;
  double t2599;
  double t2602;
  double t2603;
  double t2611;
  double t2612;
  double t2613;
  double t2614;
  double t2615;
  double t2617;
  double t2621;
  double t2622;
  double t2623;
  double t2624;
  double t2564;
  double t2570;
  double t2585;
  double t2586;
  double t2604;
  double t2605;
  double t2606;
  double t2607;
  double t2608;
  double t2609;
  double t2610;
  double t2635;
  double t2636;
  double t2637;
  double t2616;
  double t2618;
  double t2619;
  double t2620;
  double t2625;
  double t2626;
  double t2627;
  double t2628;
  double t2629;
  double t2630;
  double t2631;
  double t2646;
  double t2647;
  double t2648;
  t2524 = Sin(var1[2]);
  t2533 = Cos(var1[3]);
  t2541 = -1.*t2533*t2524;
  t2545 = Cos(var1[2]);
  t2554 = Sin(var1[3]);
  t2555 = -1.*t2545*t2554;
  t2556 = t2541 + t2555;
  t2566 = Cos(var1[4]);
  t2587 = -1.*t2545*t2533;
  t2599 = t2524*t2554;
  t2602 = t2587 + t2599;
  t2603 = Sin(var1[4]);
  t2611 = Cos(var1[5]);
  t2612 = -1.*t2611*t2524;
  t2613 = Sin(var1[5]);
  t2614 = -1.*t2545*t2613;
  t2615 = t2612 + t2614;
  t2617 = Cos(var1[6]);
  t2621 = -1.*t2545*t2611;
  t2622 = t2524*t2613;
  t2623 = t2621 + t2622;
  t2624 = Sin(var1[6]);
  t2564 = -1.0791000000000002*t2556;
  t2570 = -1.*t2566;
  t2585 = 1. + t2570;
  t2586 = 0.4*t2585*t2556;
  t2604 = -0.4*t2602*t2603;
  t2605 = t2566*t2556;
  t2606 = t2602*t2603;
  t2607 = t2605 + t2606;
  t2608 = 0.64*t2607;
  t2609 = t2586 + t2604 + t2608;
  t2610 = -9.81*t2609;
  t2635 = t2545*t2533;
  t2636 = -1.*t2524*t2554;
  t2637 = t2635 + t2636;
  t2616 = -1.0791000000000002*t2615;
  t2618 = -1.*t2617;
  t2619 = 1. + t2618;
  t2620 = 0.4*t2619*t2615;
  t2625 = -0.4*t2623*t2624;
  t2626 = t2617*t2615;
  t2627 = t2623*t2624;
  t2628 = t2626 + t2627;
  t2629 = 0.64*t2628;
  t2630 = t2620 + t2625 + t2629;
  t2631 = -9.81*t2630;
  t2646 = t2545*t2611;
  t2647 = -1.*t2524*t2613;
  t2648 = t2646 + t2647;
  p_output1[0]=0;
  p_output1[1]=-313.92;
  p_output1[2]=65.9232*t2524 + t2564 + t2610 + t2616 + t2631;
  p_output1[3]=t2564 + t2610;
  p_output1[4]=-9.81*(-0.4*t2556*t2566 + 0.4*t2603*t2637 + 0.64*(t2605 - 1.*t2603*t2637));
  p_output1[5]=t2616 + t2631;
  p_output1[6]=-9.81*(-0.4*t2615*t2617 + 0.4*t2624*t2648 + 0.64*(t2626 - 1.*t2624*t2648));
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

#include "Ge_vec_five_link_walker.hh"

namespace SymFunction
{

void Ge_vec_five_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
