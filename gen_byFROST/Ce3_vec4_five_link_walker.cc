/*
 * Automatically Generated from Mathematica.
 * Tue 24 Oct 2023 17:04:54 GMT-04:00
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
  double t2197;
  double t2251;
  double t2065;
  double t2159;
  double t1872;
  double t2151;
  double t2237;
  double t2241;
  double t2242;
  double t2249;
  double t2250;
  double t2281;
  double t2282;
  double t2291;
  double t2310;
  double t2312;
  double t2316;
  double t2317;
  double t2318;
  double t2322;
  double t2329;
  double t2252;
  double t2260;
  double t2263;
  double t2275;
  double t2278;
  double t2279;
  double t2331;
  double t2336;
  double t2337;
  double t2345;
  double t2346;
  double t2330;
  double t2332;
  double t2347;
  double t2084;
  double t2164;
  double t2171;
  double t2190;
  double t2296;
  double t2362;
  double t2363;
  double t2364;
  double t2422;
  double t2423;
  double t2424;
  double t2401;
  double t2402;
  double t2404;
  double t2409;
  double t2410;
  double t2414;
  double t2446;
  double t2450;
  double t2451;
  double t2453;
  double t2454;
  double t2455;
  double t2280;
  double t2308;
  double t2361;
  double t2365;
  double t2366;
  double t2367;
  double t2370;
  double t2371;
  double t2372;
  double t2373;
  double t2374;
  double t2383;
  double t2421;
  double t2425;
  double t2479;
  double t2480;
  double t2437;
  double t2482;
  double t2483;
  double t2439;
  t2197 = Cos(var1[4]);
  t2251 = Sin(var1[4]);
  t2065 = Sin(var1[2]);
  t2159 = Sin(var1[3]);
  t1872 = Cos(var1[3]);
  t2151 = Cos(var1[2]);
  t2237 = -1.*t2197;
  t2241 = 1. + t2237;
  t2242 = 0.4*t2241;
  t2249 = 0.64*t2197;
  t2250 = t2242 + t2249;
  t2281 = t1872*t2197;
  t2282 = -1.*t2159*t2251;
  t2291 = t2281 + t2282;
  t2310 = t2250*t2197;
  t2312 = Power(t2251,2);
  t2316 = 0.24*t2312;
  t2317 = t2310 + t2316;
  t2318 = -1.*t2197*t2159;
  t2322 = -1.*t1872*t2251;
  t2329 = t2318 + t2322;
  t2252 = t2250*t2251;
  t2260 = -0.24*t2197*t2251;
  t2263 = t2252 + t2260;
  t2275 = t2197*t2159;
  t2278 = t1872*t2251;
  t2279 = t2275 + t2278;
  t2331 = -1.*t2065*t2291;
  t2336 = -1.*t2151*t1872;
  t2337 = t2065*t2159;
  t2345 = t2336 + t2337;
  t2346 = -0.11*t2345;
  t2330 = t2151*t2329;
  t2332 = t2330 + t2331;
  t2347 = -1.*t2065*t2329;
  t2084 = -1.*t1872*t2065;
  t2164 = -1.*t2151*t2159;
  t2171 = t2084 + t2164;
  t2190 = -0.11*t2171;
  t2296 = t2151*t2291;
  t2362 = -1.*t1872*t2197;
  t2363 = t2159*t2251;
  t2364 = t2362 + t2363;
  t2422 = t1872*t2250;
  t2423 = -0.24*t2159*t2251;
  t2424 = t2422 + t2423;
  t2401 = -1.*t2250*t2159;
  t2402 = -0.24*t1872*t2251;
  t2404 = t2401 + t2402;
  t2409 = t2250*t2159;
  t2410 = 0.24*t1872*t2251;
  t2414 = t2409 + t2410;
  t2446 = -1.*t2250*t2251;
  t2450 = 0.24*t2197*t2251;
  t2451 = t2446 + t2450;
  t2453 = Power(t2197,2);
  t2454 = -0.24*t2453;
  t2455 = t2310 + t2454;
  t2280 = -1.*t2065*t2279;
  t2308 = t2280 + t2296;
  t2361 = -1.*t2263*t2332;
  t2365 = t2151*t2364;
  t2366 = t2347 + t2365;
  t2367 = -1.*t2317*t2366;
  t2370 = t2065*t2329;
  t2371 = t2370 + t2296;
  t2372 = -1.*t2263*t2371;
  t2373 = t2065*t2364;
  t2374 = t2330 + t2373;
  t2383 = -1.*t2317*t2374;
  t2421 = t2414*t2291;
  t2425 = t2329*t2424;
  t2479 = -0.24*t2197*t2159;
  t2480 = t2479 + t2402;
  t2437 = -1.*t2329*t2414;
  t2482 = 0.24*t1872*t2197;
  t2483 = t2482 + t2423;
  t2439 = -1.*t2424*t2364;
  p_output1[0]=0;
  p_output1[1]=0;
  p_output1[2]=(-0.5*(t2190 - 1.*t2263*t2308 - 1.*t2317*t2332)*var2[0] - 0.5*(-1.*t2263*(-1.*t2151*t2279 + t2331) + t2346 - 1.*t2317*(-1.*t2151*t2291 + t2347))*var2[1])*var2[3];
  p_output1[3]=(-0.5*(t2190 + t2372 + t2383)*var2[0] - 0.5*(t2346 + t2361 + t2367)*var2[1] - 0.5*(-1.*t2317*(t2291*t2404 + t2421 + t2279*t2424 + t2425) - 1.*t2263*(-1.*t2329*t2404 - 1.*t2291*t2424 + t2437 + t2439))*var2[2])*var2[3];
  p_output1[4]=var2[3]*(-0.5*(t2372 + t2383 - 1.*t2371*t2451 - 1.*(t2151*t2279 + t2065*t2291)*t2455)*var2[0] - 0.5*(t2361 + t2367 - 1.*t2332*t2451 - 1.*t2308*t2455)*var2[1] - 0.5*(-1.*(t2279*t2414 + t2291*t2424)*t2451 - 1.*(-1.*t2291*t2414 - 1.*t2329*t2424)*t2455 - 1.*t2317*(t2421 + t2425 + t2291*t2480 + t2279*t2483) - 1.*t2263*(t2437 + t2439 - 1.*t2329*t2480 - 1.*t2291*t2483))*var2[2] - 0.5*(-2.*t2317*t2451 - 2.*t2263*t2455)*var2[3] + 0.12*t2451*var2[4]);
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

#include "Ce3_vec4_five_link_walker.hh"

namespace SymFunction
{

void Ce3_vec4_five_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
