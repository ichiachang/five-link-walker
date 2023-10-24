/*
 * Automatically Generated from Mathematica.
 * Tue 24 Oct 2023 17:04:40 GMT-04:00
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
  double t1006;
  double t994;
  double t995;
  double t1007;
  double t983;
  double t996;
  double t1011;
  double t1015;
  double t1016;
  double t1017;
  double t1029;
  double t1030;
  double t1031;
  double t1032;
  double t1033;
  double t1018;
  double t1019;
  double t1020;
  double t1042;
  double t1046;
  double t1047;
  double t1063;
  double t1065;
  double t1067;
  double t1068;
  double t1069;
  double t1078;
  double t1081;
  double t1084;
  double t1071;
  double t1074;
  double t1075;
  double t1076;
  double t1077;
  double t1085;
  double t1086;
  double t1087;
  double t1088;
  t1006 = Cos(var1[5]);
  t994 = Cos(var1[6]);
  t995 = Sin(var1[5]);
  t1007 = Sin(var1[6]);
  t983 = Cos(var1[2]);
  t996 = -1.*t994*t995;
  t1011 = -1.*t1006*t1007;
  t1015 = t996 + t1011;
  t1016 = t983*t1015;
  t1017 = Sin(var1[2]);
  t1029 = -1.*t1006*t994;
  t1030 = t995*t1007;
  t1031 = t1029 + t1030;
  t1032 = t1017*t1031;
  t1033 = t1016 + t1032;
  t1018 = t1006*t994;
  t1019 = -1.*t995*t1007;
  t1020 = t1018 + t1019;
  t1042 = -1.*t1017*t1015;
  t1046 = t983*t1031;
  t1047 = t1042 + t1046;
  t1063 = -1.*t994;
  t1065 = 1. + t1063;
  t1067 = 0.4*t1065;
  t1068 = 0.64*t994;
  t1069 = t1067 + t1068;
  t1078 = t1006*t1069;
  t1081 = -0.24*t995*t1007;
  t1084 = t1078 + t1081;
  t1071 = -0.24*t1006*t1007;
  t1074 = t1069*t995;
  t1075 = 0.24*t1006*t1007;
  t1076 = t1074 + t1075;
  t1077 = t1076*t1020;
  t1085 = t1015*t1084;
  t1086 = t994*t995;
  t1087 = t1006*t1007;
  t1088 = t1086 + t1087;
  p_output1[0]=var2[6]*(-0.12*(t1016 - 1.*t1017*t1020)*var2[2] - 0.12*t1033*var2[5] - 0.12*t1033*var2[6]);
  p_output1[1]=var2[6]*(-0.12*(t1042 - 1.*t1020*t983)*var2[2] - 0.12*t1047*var2[5] - 0.12*t1047*var2[6]);
  p_output1[2]=var2[6]*(-0.12*(t1077 + t1085 + t1084*t1088 + t1020*(t1071 - 1.*t1069*t995))*var2[5] - 0.12*(t1077 + t1085 + t1088*(t1081 + 0.24*t1006*t994) + t1020*(t1071 - 0.24*t994*t995))*var2[6]);
  p_output1[3]=0;
  p_output1[4]=0;
  p_output1[5]=-0.12*(-1.*t1007*t1069 + 0.24*t1007*t994)*Power(var2[6],2);
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

#include "Ce1_vec7_five_link_walker.hh"

namespace SymFunction
{

void Ce1_vec7_five_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
