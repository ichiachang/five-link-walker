/*
 * Automatically Generated from Mathematica.
 * Tue 24 Oct 2023 17:04:36 GMT-04:00
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
  double t720;
  double t787;
  double t691;
  double t705;
  double t640;
  double t693;
  double t725;
  double t728;
  double t778;
  double t785;
  double t786;
  double t830;
  double t831;
  double t832;
  double t692;
  double t709;
  double t710;
  double t716;
  double t799;
  double t814;
  double t815;
  double t849;
  double t850;
  double t851;
  double t838;
  double t845;
  double t846;
  double t847;
  double t848;
  double t852;
  double t862;
  double t863;
  double t865;
  double t816;
  double t824;
  double t828;
  double t866;
  double t867;
  double t868;
  double t869;
  double t871;
  double t875;
  double t853;
  double t893;
  double t894;
  double t895;
  double t896;
  double t854;
  double t897;
  double t882;
  double t883;
  double t884;
  double t829;
  double t843;
  double t913;
  double t878;
  double t879;
  double t880;
  double t914;
  double t915;
  double t916;
  double t940;
  double t941;
  double t942;
  double t925;
  double t926;
  double t927;
  double t929;
  double t930;
  double t935;
  double t939;
  double t943;
  double t963;
  double t964;
  double t948;
  double t966;
  double t967;
  double t950;
  t720 = Cos(var1[4]);
  t787 = Sin(var1[4]);
  t691 = Sin(var1[2]);
  t705 = Sin(var1[3]);
  t640 = Cos(var1[3]);
  t693 = Cos(var1[2]);
  t725 = -1.*t720;
  t728 = 1. + t725;
  t778 = 0.4*t728;
  t785 = 0.64*t720;
  t786 = t778 + t785;
  t830 = t640*t720;
  t831 = -1.*t705*t787;
  t832 = t830 + t831;
  t692 = -1.*t640*t691;
  t709 = -1.*t693*t705;
  t710 = t692 + t709;
  t716 = 0.11*t710;
  t799 = t786*t787;
  t814 = -0.24*t720*t787;
  t815 = t799 + t814;
  t849 = -1.*t720*t705;
  t850 = -1.*t640*t787;
  t851 = t849 + t850;
  t838 = t693*t832;
  t845 = t786*t720;
  t846 = Power(t787,2);
  t847 = 0.24*t846;
  t848 = t845 + t847;
  t852 = t693*t851;
  t862 = t691*t851;
  t863 = t862 + t838;
  t865 = t815*t863;
  t816 = t720*t705;
  t824 = t640*t787;
  t828 = t816 + t824;
  t866 = -1.*t640*t720;
  t867 = t705*t787;
  t868 = t866 + t867;
  t869 = t691*t868;
  t871 = t852 + t869;
  t875 = t848*t871;
  t853 = -1.*t691*t832;
  t893 = -1.*t693*t640;
  t894 = t691*t705;
  t895 = t893 + t894;
  t896 = 0.11*t895;
  t854 = t852 + t853;
  t897 = -1.*t691*t851;
  t882 = Power(t720,2);
  t883 = -0.24*t882;
  t884 = t845 + t883;
  t829 = -1.*t691*t828;
  t843 = t829 + t838;
  t913 = t815*t854;
  t878 = -1.*t786*t787;
  t879 = 0.24*t720*t787;
  t880 = t878 + t879;
  t914 = t693*t868;
  t915 = t897 + t914;
  t916 = t848*t915;
  t940 = t640*t786;
  t941 = -0.24*t705*t787;
  t942 = t940 + t941;
  t925 = -1.*t786*t705;
  t926 = -0.24*t640*t787;
  t927 = t925 + t926;
  t929 = t786*t705;
  t930 = 0.24*t640*t787;
  t935 = t929 + t930;
  t939 = t935*t832;
  t943 = t851*t942;
  t963 = -0.24*t720*t705;
  t964 = t963 + t926;
  t948 = -1.*t851*t935;
  t966 = 0.24*t640*t720;
  t967 = t966 + t941;
  t950 = -1.*t942*t868;
  p_output1[0]=var2[3]*(-0.5*(t716 + t815*t843 + t848*t854)*var2[2] - 0.5*(t716 + t865 + t875)*var2[3] - 0.5*(t865 + t875 + t863*t880 + (t693*t828 + t691*t832)*t884)*var2[4]);
  p_output1[1]=var2[3]*(-0.5*(t815*(-1.*t693*t828 + t853) + t896 + t848*(-1.*t693*t832 + t897))*var2[2] - 0.5*(t896 + t913 + t916)*var2[3] - 0.5*(t854*t880 + t843*t884 + t913 + t916)*var2[4]);
  p_output1[2]=var2[3]*(-0.5*(t848*(t832*t927 + t939 + t828*t942 + t943) + t815*(-1.*t851*t927 - 1.*t832*t942 + t948 + t950))*var2[3] - 0.5*(t880*(t828*t935 + t832*t942) + t884*(-1.*t832*t935 - 1.*t851*t942) + t848*(t939 + t943 + t832*t964 + t828*t967) + t815*(t948 + t950 - 1.*t851*t964 - 1.*t832*t967))*var2[4]);
  p_output1[3]=-0.5*(2.*t848*t880 + 2.*t815*t884)*var2[3]*var2[4];
  p_output1[4]=-0.12*t880*var2[3]*var2[4];
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

#include "Ce1_vec4_five_link_walker.hh"

namespace SymFunction
{

void Ce1_vec4_five_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
