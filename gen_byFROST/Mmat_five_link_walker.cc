/*
 * Automatically Generated from Mathematica.
 * Tue 24 Oct 2023 17:04:29 GMT-04:00
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
static void output1(double *p_output1,const double *var1)
{
  double t117;
  double t112;
  double t120;
  double t122;
  double t130;
  double t132;
  double t136;
  double t137;
  double t138;
  double t149;
  double t151;
  double t159;
  double t161;
  double t165;
  double t166;
  double t167;
  double t126;
  double t127;
  double t128;
  double t121;
  double t123;
  double t124;
  double t131;
  double t133;
  double t134;
  double t135;
  double t139;
  double t140;
  double t142;
  double t143;
  double t144;
  double t145;
  double t146;
  double t147;
  double t155;
  double t156;
  double t157;
  double t150;
  double t152;
  double t153;
  double t160;
  double t162;
  double t163;
  double t164;
  double t168;
  double t169;
  double t171;
  double t172;
  double t173;
  double t174;
  double t175;
  double t176;
  double t211;
  double t212;
  double t213;
  double t214;
  double t215;
  double t216;
  double t217;
  double t218;
  double t220;
  double t221;
  double t222;
  double t236;
  double t237;
  double t238;
  double t239;
  double t240;
  double t241;
  double t242;
  double t243;
  double t245;
  double t246;
  double t247;
  double t179;
  double t180;
  double t181;
  double t182;
  double t183;
  double t184;
  double t185;
  double t186;
  double t187;
  double t188;
  double t189;
  double t190;
  double t191;
  double t192;
  double t193;
  double t194;
  double t195;
  double t196;
  double t197;
  double t198;
  double t199;
  double t200;
  double t201;
  double t202;
  double t203;
  double t115;
  double t116;
  double t118;
  double t119;
  double t129;
  double t158;
  double t205;
  double t206;
  double t207;
  double t208;
  double t209;
  double t219;
  double t223;
  double t224;
  double t226;
  double t227;
  double t228;
  double t230;
  double t231;
  double t232;
  double t233;
  double t234;
  double t244;
  double t248;
  double t249;
  double t251;
  double t252;
  double t253;
  double t262;
  double t263;
  double t264;
  double t257;
  double t258;
  double t259;
  double t260;
  double t274;
  double t275;
  double t276;
  double t269;
  double t270;
  double t271;
  double t272;
  double t204;
  double t210;
  double t225;
  double t229;
  double t235;
  double t250;
  double t254;
  double t255;
  double t287;
  double t288;
  double t289;
  double t290;
  double t291;
  double t292;
  double t293;
  double t294;
  double t256;
  double t261;
  double t265;
  double t266;
  double t295;
  double t296;
  double t297;
  double t298;
  double t312;
  double t313;
  double t314;
  double t315;
  double t267;
  double t299;
  double t316;
  double t317;
  double t327;
  double t328;
  double t268;
  double t273;
  double t277;
  double t278;
  double t300;
  double t301;
  double t302;
  double t303;
  double t318;
  double t319;
  double t320;
  double t321;
  double t279;
  double t304;
  double t322;
  double t323;
  double t332;
  double t333;
  t117 = Sin(var1[2]);
  t112 = Cos(var1[2]);
  t120 = Cos(var1[3]);
  t122 = Sin(var1[3]);
  t130 = Cos(var1[4]);
  t132 = Sin(var1[4]);
  t136 = t120*t130;
  t137 = -1.*t122*t132;
  t138 = t136 + t137;
  t149 = Cos(var1[5]);
  t151 = Sin(var1[5]);
  t159 = Cos(var1[6]);
  t161 = Sin(var1[6]);
  t165 = t149*t159;
  t166 = -1.*t151*t161;
  t167 = t165 + t166;
  t126 = t112*t120;
  t127 = -1.*t117*t122;
  t128 = t126 + t127;
  t121 = t120*t117;
  t123 = t112*t122;
  t124 = t121 + t123;
  t131 = -1.*t130*t122;
  t133 = -1.*t120*t132;
  t134 = t131 + t133;
  t135 = t117*t134;
  t139 = t112*t138;
  t140 = t135 + t139;
  t142 = t130*t122;
  t143 = t120*t132;
  t144 = t142 + t143;
  t145 = t112*t144;
  t146 = t117*t138;
  t147 = t145 + t146;
  t155 = t112*t149;
  t156 = -1.*t117*t151;
  t157 = t155 + t156;
  t150 = t149*t117;
  t152 = t112*t151;
  t153 = t150 + t152;
  t160 = -1.*t159*t151;
  t162 = -1.*t149*t161;
  t163 = t160 + t162;
  t164 = t117*t163;
  t168 = t112*t167;
  t169 = t164 + t168;
  t171 = t159*t151;
  t172 = t149*t161;
  t173 = t171 + t172;
  t174 = t112*t173;
  t175 = t117*t167;
  t176 = t174 + t175;
  t211 = -1.*t130;
  t212 = 1. + t211;
  t213 = 0.4*t212;
  t214 = 0.64*t130;
  t215 = t213 + t214;
  t216 = t215*t122;
  t217 = 0.24*t120*t132;
  t218 = t216 + t217;
  t220 = t120*t215;
  t221 = -0.24*t122*t132;
  t222 = t220 + t221;
  t236 = -1.*t159;
  t237 = 1. + t236;
  t238 = 0.4*t237;
  t239 = 0.64*t159;
  t240 = t238 + t239;
  t241 = t240*t151;
  t242 = 0.24*t149*t161;
  t243 = t241 + t242;
  t245 = t149*t240;
  t246 = -0.24*t151*t161;
  t247 = t245 + t246;
  t179 = -1.*t120*t117;
  t180 = -1.*t112*t122;
  t181 = t179 + t180;
  t182 = t181*t128;
  t183 = t124*t128;
  t184 = t112*t134;
  t185 = -1.*t117*t138;
  t186 = t184 + t185;
  t187 = t140*t186;
  t188 = -1.*t117*t144;
  t189 = t188 + t139;
  t190 = t189*t147;
  t191 = -1.*t149*t117;
  t192 = -1.*t112*t151;
  t193 = t191 + t192;
  t194 = t193*t157;
  t195 = t153*t157;
  t196 = t112*t163;
  t197 = -1.*t117*t167;
  t198 = t196 + t197;
  t199 = t169*t198;
  t200 = -1.*t117*t173;
  t201 = t200 + t168;
  t202 = t201*t176;
  t203 = t182 + t183 + t187 + t190 + t194 + t195 + t199 + t202;
  t115 = Power(t112,2);
  t116 = 28.*t115;
  t118 = Power(t117,2);
  t119 = 28.*t118;
  t129 = Power(t128,2);
  t158 = Power(t157,2);
  t205 = Power(t120,2);
  t206 = 0.11*t205;
  t207 = Power(t122,2);
  t208 = 0.11*t207;
  t209 = t206 + t208;
  t219 = -1.*t218*t138;
  t223 = -1.*t134*t222;
  t224 = t219 + t223;
  t226 = t218*t144;
  t227 = t138*t222;
  t228 = t226 + t227;
  t230 = Power(t149,2);
  t231 = 0.11*t230;
  t232 = Power(t151,2);
  t233 = 0.11*t232;
  t234 = t231 + t233;
  t244 = -1.*t243*t167;
  t248 = -1.*t163*t247;
  t249 = t244 + t248;
  t251 = t243*t173;
  t252 = t167*t247;
  t253 = t251 + t252;
  t262 = t215*t132;
  t263 = -0.24*t130*t132;
  t264 = t262 + t263;
  t257 = t215*t130;
  t258 = Power(t132,2);
  t259 = 0.24*t258;
  t260 = t257 + t259;
  t274 = t240*t161;
  t275 = -0.24*t159*t161;
  t276 = t274 + t275;
  t269 = t240*t159;
  t270 = Power(t161,2);
  t271 = 0.24*t270;
  t272 = t269 + t271;
  t204 = 6.72*t112;
  t210 = t128*t209;
  t225 = t147*t224;
  t229 = t140*t228;
  t235 = t157*t234;
  t250 = t176*t249;
  t254 = t169*t253;
  t255 = t204 + t210 + t225 + t229 + t235 + t250 + t254;
  t287 = -6.72*t117;
  t288 = t181*t209;
  t289 = t189*t224;
  t290 = t186*t228;
  t291 = t193*t234;
  t292 = t201*t249;
  t293 = t198*t253;
  t294 = t287 + t288 + t289 + t290 + t291 + t292 + t293;
  t256 = 0.11*t128;
  t261 = t260*t140;
  t265 = t264*t147;
  t266 = t256 + t261 + t265;
  t295 = 0.11*t181;
  t296 = t264*t189;
  t297 = t260*t186;
  t298 = t295 + t296 + t297;
  t312 = 0.11*t209;
  t313 = t264*t224;
  t314 = t260*t228;
  t315 = 0.67 + t312 + t313 + t314;
  t267 = 0.24*t140;
  t299 = 0.24*t186;
  t316 = 0.24*t228;
  t317 = 0.2 + t316;
  t327 = 0.24*t260;
  t328 = 0.2 + t327;
  t268 = 0.11*t157;
  t273 = t272*t169;
  t277 = t276*t176;
  t278 = t268 + t273 + t277;
  t300 = 0.11*t193;
  t301 = t276*t201;
  t302 = t272*t198;
  t303 = t300 + t301 + t302;
  t318 = 0.11*t234;
  t319 = t276*t249;
  t320 = t272*t253;
  t321 = 0.67 + t318 + t319 + t320;
  t279 = 0.24*t169;
  t304 = 0.24*t198;
  t322 = 0.24*t253;
  t323 = 0.2 + t322;
  t332 = 0.24*t272;
  t333 = 0.2 + t332;
  p_output1[0]=t116 + t119 + Power(t124,2) + t129 + Power(t140,2) + Power(t147,2) + Power(t153,2) + t158 + Power(t169,2) + Power(t176,2);
  p_output1[1]=t203;
  p_output1[2]=t255;
  p_output1[3]=t266;
  p_output1[4]=t267;
  p_output1[5]=t278;
  p_output1[6]=t279;
  p_output1[7]=t203;
  p_output1[8]=t116 + t119 + t129 + t158 + Power(t181,2) + Power(t186,2) + Power(t189,2) + Power(t193,2) + Power(t198,2) + Power(t201,2);
  p_output1[9]=t294;
  p_output1[10]=t298;
  p_output1[11]=t299;
  p_output1[12]=t303;
  p_output1[13]=t304;
  p_output1[14]=t255;
  p_output1[15]=t294;
  p_output1[16]=4.2828 + Power(t209,2) + Power(t224,2) + Power(t228,2) + Power(t234,2) + Power(t249,2) + Power(t253,2);
  p_output1[17]=t315;
  p_output1[18]=t317;
  p_output1[19]=t321;
  p_output1[20]=t323;
  p_output1[21]=t266;
  p_output1[22]=t298;
  p_output1[23]=t315;
  p_output1[24]=1.5121 + Power(t260,2) + Power(t264,2);
  p_output1[25]=t328;
  p_output1[26]=0;
  p_output1[27]=0;
  p_output1[28]=t267;
  p_output1[29]=t299;
  p_output1[30]=t317;
  p_output1[31]=t328;
  p_output1[32]=1.0876;
  p_output1[33]=0;
  p_output1[34]=0;
  p_output1[35]=t278;
  p_output1[36]=t303;
  p_output1[37]=t321;
  p_output1[38]=0;
  p_output1[39]=0;
  p_output1[40]=1.5121 + Power(t272,2) + Power(t276,2);
  p_output1[41]=t333;
  p_output1[42]=t279;
  p_output1[43]=t304;
  p_output1[44]=t323;
  p_output1[45]=0;
  p_output1[46]=0;
  p_output1[47]=t333;
  p_output1[48]=1.0876;
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

  double *var1;
  double *p_output1;

  /*  Check for proper number of arguments.  */ 
  if( nrhs != 1)
    {
      mexErrMsgIdAndTxt("MATLAB:MShaped:invalidNumInputs", "One input(s) required (var1).");
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

  /*  Assign pointers to each input.  */
  var1 = mxGetPr(prhs[0]);
   


   
  /*  Create matrices for return arguments.  */
  plhs[0] = mxCreateDoubleMatrix((mwSize) 7, (mwSize) 7, mxREAL);
  p_output1 = mxGetPr(plhs[0]);


  /* Call the calculation subroutine. */
  output1(p_output1,var1);


}

#else // MATLAB_MEX_FILE

#include "Mmat_five_link_walker.hh"

namespace SymFunction
{

void Mmat_five_link_walker_raw(double *p_output1, const double *var1)
{
  // Call Subroutines
  output1(p_output1, var1);

}

}

#endif // MATLAB_MEX_FILE
