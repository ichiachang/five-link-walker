/*
 * Automatically Generated from Mathematica.
 * Tue 24 Oct 2023 17:04:39 GMT-04:00
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
  double t921;
  double t954;
  double t892;
  double t912;
  double t891;
  double t911;
  double t922;
  double t928;
  double t944;
  double t945;
  double t946;
  double t978;
  double t979;
  double t980;
  double t898;
  double t917;
  double t918;
  double t920;
  double t956;
  double t957;
  double t973;
  double t988;
  double t989;
  double t990;
  double t981;
  double t984;
  double t985;
  double t986;
  double t987;
  double t991;
  double t997;
  double t998;
  double t999;
  double t974;
  double t975;
  double t976;
  double t1000;
  double t1001;
  double t1002;
  double t1003;
  double t1004;
  double t1005;
  double t992;
  double t1023;
  double t1024;
  double t1025;
  double t1026;
  double t993;
  double t1027;
  double t1012;
  double t1013;
  double t1014;
  double t977;
  double t982;
  double t1036;
  double t1008;
  double t1009;
  double t1010;
  double t1037;
  double t1038;
  double t1039;
  double t1056;
  double t1057;
  double t1058;
  double t1048;
  double t1049;
  double t1050;
  double t1052;
  double t1053;
  double t1054;
  double t1055;
  double t1059;
  double t1079;
  double t1080;
  double t1064;
  double t1082;
  double t1083;
  double t1066;
  t921 = Cos(var1[6]);
  t954 = Sin(var1[6]);
  t892 = Sin(var1[2]);
  t912 = Sin(var1[5]);
  t891 = Cos(var1[5]);
  t911 = Cos(var1[2]);
  t922 = -1.*t921;
  t928 = 1. + t922;
  t944 = 0.4*t928;
  t945 = 0.64*t921;
  t946 = t944 + t945;
  t978 = t891*t921;
  t979 = -1.*t912*t954;
  t980 = t978 + t979;
  t898 = -1.*t891*t892;
  t917 = -1.*t911*t912;
  t918 = t898 + t917;
  t920 = 0.11*t918;
  t956 = t946*t954;
  t957 = -0.24*t921*t954;
  t973 = t956 + t957;
  t988 = -1.*t921*t912;
  t989 = -1.*t891*t954;
  t990 = t988 + t989;
  t981 = t911*t980;
  t984 = t946*t921;
  t985 = Power(t954,2);
  t986 = 0.24*t985;
  t987 = t984 + t986;
  t991 = t911*t990;
  t997 = t892*t990;
  t998 = t997 + t981;
  t999 = t973*t998;
  t974 = t921*t912;
  t975 = t891*t954;
  t976 = t974 + t975;
  t1000 = -1.*t891*t921;
  t1001 = t912*t954;
  t1002 = t1000 + t1001;
  t1003 = t892*t1002;
  t1004 = t991 + t1003;
  t1005 = t987*t1004;
  t992 = -1.*t892*t980;
  t1023 = -1.*t911*t891;
  t1024 = t892*t912;
  t1025 = t1023 + t1024;
  t1026 = 0.11*t1025;
  t993 = t991 + t992;
  t1027 = -1.*t892*t990;
  t1012 = Power(t921,2);
  t1013 = -0.24*t1012;
  t1014 = t984 + t1013;
  t977 = -1.*t892*t976;
  t982 = t977 + t981;
  t1036 = t973*t993;
  t1008 = -1.*t946*t954;
  t1009 = 0.24*t921*t954;
  t1010 = t1008 + t1009;
  t1037 = t911*t1002;
  t1038 = t1027 + t1037;
  t1039 = t987*t1038;
  t1056 = t891*t946;
  t1057 = -0.24*t912*t954;
  t1058 = t1056 + t1057;
  t1048 = -1.*t946*t912;
  t1049 = -0.24*t891*t954;
  t1050 = t1048 + t1049;
  t1052 = t946*t912;
  t1053 = 0.24*t891*t954;
  t1054 = t1052 + t1053;
  t1055 = t1054*t980;
  t1059 = t990*t1058;
  t1079 = -0.24*t921*t912;
  t1080 = t1079 + t1049;
  t1064 = -1.*t990*t1054;
  t1082 = 0.24*t891*t921;
  t1083 = t1082 + t1057;
  t1066 = -1.*t1058*t1002;
  p_output1[0]=var2[5]*(-0.5*(t920 + t973*t982 + t987*t993)*var2[2] - 0.5*(t1005 + t920 + t999)*var2[5] - 0.5*(t1005 + t1014*(t911*t976 + t892*t980) + t1010*t998 + t999)*var2[6]);
  p_output1[1]=var2[5]*(-0.5*(t1026 + (t1027 - 1.*t911*t980)*t987 + t973*(-1.*t911*t976 + t992))*var2[2] - 0.5*(t1026 + t1036 + t1039)*var2[5] - 0.5*(t1036 + t1039 + t1014*t982 + t1010*t993)*var2[6]);
  p_output1[2]=var2[5]*(-0.5*((t1055 + t1059 + t1058*t976 + t1050*t980)*t987 + t973*(t1064 + t1066 - 1.*t1058*t980 - 1.*t1050*t990))*var2[5] - 0.5*(t1010*(t1054*t976 + t1058*t980) + (t1055 + t1059 + t1083*t976 + t1080*t980)*t987 + t1014*(-1.*t1054*t980 - 1.*t1058*t990) + t973*(t1064 + t1066 - 1.*t1083*t980 - 1.*t1080*t990))*var2[6]);
  p_output1[3]=0;
  p_output1[4]=0;
  p_output1[5]=-0.5*(2.*t1014*t973 + 2.*t1010*t987)*var2[5]*var2[6];
  p_output1[6]=-0.12*t1010*var2[5]*var2[6];
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

#include "Ce1_vec6_five_link_walker.hh"

namespace SymFunction
{

void Ce1_vec6_five_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
