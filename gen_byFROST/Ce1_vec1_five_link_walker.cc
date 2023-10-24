/*
 * Automatically Generated from Mathematica.
 * Tue 24 Oct 2023 17:04:30 GMT-04:00
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
  double t126;
  double t123;
  double t124;
  double t129;
  double t141;
  double t121;
  double t142;
  double t143;
  double t145;
  double t125;
  double t135;
  double t136;
  double t139;
  double t146;
  double t147;
  double t183;
  double t184;
  double t187;
  double t148;
  double t150;
  double t152;
  double t153;
  double t154;
  double t155;
  double t158;
  double t164;
  double t165;
  double t168;
  double t170;
  double t171;
  double t172;
  double t174;
  double t177;
  double t178;
  double t182;
  double t189;
  double t190;
  double t194;
  double t195;
  double t196;
  double t225;
  double t227;
  double t240;
  double t241;
  double t245;
  double t257;
  double t262;
  double t261;
  double t264;
  double t265;
  double t269;
  double t273;
  double t274;
  double t276;
  double t283;
  double t284;
  double t285;
  double t266;
  double t277;
  double t289;
  double t290;
  double t291;
  double t278;
  double t226;
  double t229;
  double t235;
  double t250;
  double t251;
  double t252;
  double t254;
  double t255;
  double t296;
  double t297;
  double t298;
  double t301;
  double t302;
  double t303;
  double t305;
  double t215;
  double t216;
  double t202;
  double t203;
  double t324;
  double t325;
  double t326;
  double t329;
  double t330;
  double t331;
  double t334;
  double t337;
  double t338;
  double t339;
  double t340;
  double t341;
  double t342;
  double t343;
  double t286;
  double t288;
  double t280;
  double t281;
  double t361;
  double t354;
  double t355;
  double t356;
  double t357;
  double t358;
  double t359;
  double t360;
  double t372;
  double t373;
  double t374;
  double t375;
  double t376;
  double t377;
  double t390;
  double t391;
  double t392;
  double t393;
  double t394;
  double t395;
  double t396;
  double t397;
  double t399;
  double t400;
  double t401;
  double t405;
  double t406;
  double t407;
  double t398;
  double t402;
  double t403;
  double t404;
  double t409;
  double t410;
  double t414;
  double t415;
  double t416;
  double t417;
  double t426;
  double t427;
  double t419;
  double t429;
  double t430;
  double t421;
  double t384;
  double t385;
  double t386;
  double t387;
  double t388;
  double t389;
  double t449;
  double t450;
  double t451;
  double t452;
  double t453;
  double t454;
  double t455;
  double t456;
  double t458;
  double t459;
  double t460;
  double t443;
  double t444;
  double t445;
  double t446;
  double t447;
  double t448;
  double t457;
  double t461;
  double t462;
  double t464;
  double t465;
  double t466;
  double t471;
  double t472;
  double t473;
  double t470;
  double t475;
  double t476;
  double t480;
  double t489;
  double t490;
  double t482;
  double t492;
  double t493;
  double t484;
  double t505;
  double t506;
  double t507;
  double t508;
  double t510;
  double t511;
  double t512;
  double t513;
  double t517;
  double t518;
  double t538;
  double t539;
  double t540;
  double t541;
  double t543;
  double t544;
  double t545;
  double t546;
  double t550;
  double t551;
  t126 = Cos(var1[3]);
  t123 = Cos(var1[4]);
  t124 = Sin(var1[3]);
  t129 = Sin(var1[4]);
  t141 = Cos(var1[2]);
  t121 = Sin(var1[2]);
  t142 = t126*t123;
  t143 = -1.*t124*t129;
  t145 = t142 + t143;
  t125 = -1.*t123*t124;
  t135 = -1.*t126*t129;
  t136 = t125 + t135;
  t139 = t121*t136;
  t146 = t141*t145;
  t147 = t139 + t146;
  t183 = t141*t126;
  t184 = -1.*t121*t124;
  t187 = t183 + t184;
  t148 = t123*t124;
  t150 = t126*t129;
  t152 = t148 + t150;
  t153 = t141*t152;
  t154 = t121*t145;
  t155 = t153 + t154;
  t158 = 2.*t147*t155;
  t164 = t141*t136;
  t165 = -1.*t126*t123;
  t168 = t124*t129;
  t170 = t165 + t168;
  t171 = t121*t170;
  t172 = t164 + t171;
  t174 = 2.*t147*t172;
  t177 = -1.*t126*t121;
  t178 = -1.*t141*t124;
  t182 = t177 + t178;
  t189 = 2.*t182*t187;
  t190 = t126*t121;
  t194 = t141*t124;
  t195 = t190 + t194;
  t196 = 2.*t195*t187;
  t225 = Cos(var1[5]);
  t227 = Sin(var1[5]);
  t240 = t141*t225;
  t241 = -1.*t121*t227;
  t245 = t240 + t241;
  t257 = Cos(var1[6]);
  t262 = Sin(var1[6]);
  t261 = -1.*t257*t227;
  t264 = -1.*t225*t262;
  t265 = t261 + t264;
  t269 = t225*t257;
  t273 = -1.*t227*t262;
  t274 = t269 + t273;
  t276 = t141*t274;
  t283 = t257*t227;
  t284 = t225*t262;
  t285 = t283 + t284;
  t266 = t121*t265;
  t277 = t266 + t276;
  t289 = t141*t285;
  t290 = t121*t274;
  t291 = t289 + t290;
  t278 = t141*t265;
  t226 = -1.*t225*t121;
  t229 = -1.*t141*t227;
  t235 = t226 + t229;
  t250 = 2.*t235*t245;
  t251 = t225*t121;
  t252 = t141*t227;
  t254 = t251 + t252;
  t255 = 2.*t254*t245;
  t296 = 2.*t277*t291;
  t297 = -1.*t225*t257;
  t298 = t227*t262;
  t301 = t297 + t298;
  t302 = t121*t301;
  t303 = t278 + t302;
  t305 = 2.*t277*t303;
  t215 = -1.*t121*t152;
  t216 = t215 + t146;
  t202 = -1.*t121*t145;
  t203 = t164 + t202;
  t324 = t147*t216;
  t325 = t203*t155;
  t326 = -1.*t121*t136;
  t329 = t141*t170;
  t330 = t326 + t329;
  t331 = t147*t330;
  t334 = t203*t172;
  t337 = Power(t182,2);
  t338 = t182*t195;
  t339 = Power(t187,2);
  t340 = -1.*t141*t126;
  t341 = t121*t124;
  t342 = t340 + t341;
  t343 = t187*t342;
  t286 = -1.*t121*t285;
  t288 = t286 + t276;
  t280 = -1.*t121*t274;
  t281 = t278 + t280;
  t361 = -1.*t121*t265;
  t354 = Power(t235,2);
  t355 = t235*t254;
  t356 = Power(t245,2);
  t357 = -1.*t141*t225;
  t358 = t121*t227;
  t359 = t357 + t358;
  t360 = t245*t359;
  t372 = t277*t288;
  t373 = t281*t291;
  t374 = t141*t301;
  t375 = t361 + t374;
  t376 = t277*t375;
  t377 = t281*t303;
  t390 = -1.*t123;
  t391 = 1. + t390;
  t392 = 0.4*t391;
  t393 = 0.64*t123;
  t394 = t392 + t393;
  t395 = t394*t124;
  t396 = 0.24*t126*t129;
  t397 = t395 + t396;
  t399 = t126*t394;
  t400 = -0.24*t124*t129;
  t401 = t399 + t400;
  t405 = -1.*t394*t124;
  t406 = -0.24*t126*t129;
  t407 = t405 + t406;
  t398 = -1.*t397*t145;
  t402 = -1.*t136*t401;
  t403 = t398 + t402;
  t404 = t147*t403;
  t409 = t397*t145;
  t410 = t136*t401;
  t414 = t397*t152;
  t415 = t145*t401;
  t416 = t414 + t415;
  t417 = t416*t172;
  t426 = -0.24*t123*t124;
  t427 = t426 + t406;
  t419 = -1.*t136*t397;
  t429 = 0.24*t126*t123;
  t430 = t429 + t400;
  t421 = -1.*t401*t170;
  t384 = Power(t126,2);
  t385 = 0.11*t384;
  t386 = Power(t124,2);
  t387 = 0.11*t386;
  t388 = t385 + t387;
  t389 = t182*t388;
  t449 = -1.*t257;
  t450 = 1. + t449;
  t451 = 0.4*t450;
  t452 = 0.64*t257;
  t453 = t451 + t452;
  t454 = t453*t227;
  t455 = 0.24*t225*t262;
  t456 = t454 + t455;
  t458 = t225*t453;
  t459 = -0.24*t227*t262;
  t460 = t458 + t459;
  t443 = Power(t225,2);
  t444 = 0.11*t443;
  t445 = Power(t227,2);
  t446 = 0.11*t445;
  t447 = t444 + t446;
  t448 = t235*t447;
  t457 = -1.*t456*t274;
  t461 = -1.*t265*t460;
  t462 = t457 + t461;
  t464 = t456*t285;
  t465 = t274*t460;
  t466 = t464 + t465;
  t471 = -1.*t453*t227;
  t472 = -0.24*t225*t262;
  t473 = t471 + t472;
  t470 = t277*t462;
  t475 = t456*t274;
  t476 = t265*t460;
  t480 = t466*t303;
  t489 = -0.24*t257*t227;
  t490 = t489 + t472;
  t482 = -1.*t265*t456;
  t492 = 0.24*t225*t257;
  t493 = t492 + t459;
  t484 = -1.*t460*t301;
  t505 = 0.11*t182;
  t506 = t394*t129;
  t507 = -0.24*t123*t129;
  t508 = t506 + t507;
  t510 = t394*t123;
  t511 = Power(t129,2);
  t512 = 0.24*t511;
  t513 = t510 + t512;
  t517 = t508*t147;
  t518 = t513*t172;
  t538 = 0.11*t235;
  t539 = t453*t262;
  t540 = -0.24*t257*t262;
  t541 = t539 + t540;
  t543 = t453*t257;
  t544 = Power(t262,2);
  t545 = 0.24*t544;
  t546 = t543 + t545;
  t550 = t541*t277;
  t551 = t546*t303;
  p_output1[0]=var2[0]*(-0.5*(t189 + t196 + 2.*t147*t203 + 2.*t155*t216 + t250 + t255 + 2.*t277*t281 + 2.*t288*t291)*var2[2] - 0.5*(t158 + t174 + t189 + t196)*var2[3] - 0.5*(t158 + t174)*var2[4] - 0.5*(t250 + t255 + t296 + t305)*var2[5] - 0.5*(t296 + t305)*var2[6]);
  p_output1[1]=var2[0]*(-0.5*(t155*(-1.*t141*t152 + t202) + Power(t203,2) + Power(t216,2) + Power(t281,2) + Power(t288,2) + (t280 - 1.*t141*t285)*t291 + t147*(-1.*t141*t145 + t326) + t337 + t338 + t339 + t343 + t354 + t355 + t356 + t360 + t277*(-1.*t141*t274 + t361))*var2[2] - 0.5*(t324 + t325 + t331 + t334 + t337 + t338 + t339 + t343)*var2[3] - 0.5*(t324 + t325 + t331 + t334)*var2[4] - 0.5*(t354 + t355 + t356 + t360 + t372 + t373 + t376 + t377)*var2[5] - 0.5*(t372 + t373 + t376 + t377)*var2[6]);
  p_output1[2]=var2[0]*(-0.5*(-6.72*t121 + t389 + t216*t403 + t203*t416 + t448 + t288*t462 + t281*t466)*var2[2] - 0.5*(t389 + t404 + t147*(t152*t401 + t145*t407 + t409 + t410) + t417 + t155*(-1.*t145*t401 - 1.*t136*t407 + t419 + t421))*var2[3] - 0.5*(t404 + t417 + t155*(t419 + t421 - 1.*t136*t427 - 1.*t145*t430) + t147*(t409 + t410 + t145*t427 + t152*t430))*var2[4] - 0.5*(t448 + t470 + t277*(t285*t460 + t274*t473 + t475 + t476) + t480 + t291*(-1.*t274*t460 - 1.*t265*t473 + t482 + t484))*var2[5] - 0.5*(t470 + t480 + t291*(t482 + t484 - 1.*t265*t490 - 1.*t274*t493) + t277*(t475 + t476 + t274*t490 + t285*t493))*var2[6]);
  p_output1[3]=var2[0]*(-0.5*(t505 + t216*t508 + t203*t513)*var2[2] - 0.5*(t505 + t517 + t518)*var2[3] - 0.5*(t147*(0.24*t123*t129 - 1.*t129*t394) + t155*(-0.24*Power(t123,2) + t510) + t517 + t518)*var2[4]);
  p_output1[4]=var2[0]*(-0.12*t203*var2[2] - 0.12*t172*var2[3] - 0.12*t172*var2[4]);
  p_output1[5]=var2[0]*(-0.5*(t538 + t288*t541 + t281*t546)*var2[2] - 0.5*(t538 + t550 + t551)*var2[5] - 0.5*(t277*(0.24*t257*t262 - 1.*t262*t453) + t291*(-0.24*Power(t257,2) + t543) + t550 + t551)*var2[6]);
  p_output1[6]=var2[0]*(-0.12*t281*var2[2] - 0.12*t303*var2[5] - 0.12*t303*var2[6]);
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

#include "Ce1_vec1_five_link_walker.hh"

namespace SymFunction
{

void Ce1_vec1_five_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
