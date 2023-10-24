/*
 * Automatically Generated from Mathematica.
 * Tue 24 Oct 2023 17:04:37 GMT-04:00
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
  double t876;
  double t855;
  double t856;
  double t877;
  double t844;
  double t861;
  double t881;
  double t885;
  double t886;
  double t887;
  double t899;
  double t900;
  double t901;
  double t906;
  double t910;
  double t888;
  double t889;
  double t890;
  double t919;
  double t923;
  double t924;
  double t947;
  double t949;
  double t951;
  double t952;
  double t953;
  double t962;
  double t965;
  double t968;
  double t955;
  double t958;
  double t959;
  double t960;
  double t961;
  double t969;
  double t970;
  double t971;
  double t972;
  t876 = Cos(var1[3]);
  t855 = Cos(var1[4]);
  t856 = Sin(var1[3]);
  t877 = Sin(var1[4]);
  t844 = Cos(var1[2]);
  t861 = -1.*t855*t856;
  t881 = -1.*t876*t877;
  t885 = t861 + t881;
  t886 = t844*t885;
  t887 = Sin(var1[2]);
  t899 = -1.*t876*t855;
  t900 = t856*t877;
  t901 = t899 + t900;
  t906 = t887*t901;
  t910 = t886 + t906;
  t888 = t876*t855;
  t889 = -1.*t856*t877;
  t890 = t888 + t889;
  t919 = -1.*t887*t885;
  t923 = t844*t901;
  t924 = t919 + t923;
  t947 = -1.*t855;
  t949 = 1. + t947;
  t951 = 0.4*t949;
  t952 = 0.64*t855;
  t953 = t951 + t952;
  t962 = t876*t953;
  t965 = -0.24*t856*t877;
  t968 = t962 + t965;
  t955 = -0.24*t876*t877;
  t958 = t953*t856;
  t959 = 0.24*t876*t877;
  t960 = t958 + t959;
  t961 = t960*t890;
  t969 = t885*t968;
  t970 = t855*t856;
  t971 = t876*t877;
  t972 = t970 + t971;
  p_output1[0]=var2[4]*(-0.12*(t886 - 1.*t887*t890)*var2[2] - 0.12*t910*var2[3] - 0.12*t910*var2[4]);
  p_output1[1]=var2[4]*(-0.12*(-1.*t844*t890 + t919)*var2[2] - 0.12*t924*var2[3] - 0.12*t924*var2[4]);
  p_output1[2]=var2[4]*(-0.12*(t890*(-1.*t856*t953 + t955) + t961 + t969 + t968*t972)*var2[3] - 0.12*(t890*(-0.24*t855*t856 + t955) + t961 + t969 + (0.24*t855*t876 + t965)*t972)*var2[4]);
  p_output1[3]=-0.12*(0.24*t855*t877 - 1.*t877*t953)*Power(var2[4],2);
  p_output1[4]=0;
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

#include "Ce1_vec5_five_link_walker.hh"

namespace SymFunction
{

void Ce1_vec5_five_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
