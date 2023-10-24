/*
 * Automatically Generated from Mathematica.
 * Tue 24 Oct 2023 17:04:52 GMT-04:00
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
  double t1838;
  double t1897;
  double t1956;
  double t1917;
  double t2045;
  double t2062;
  double t2066;
  double t2070;
  double t2082;
  double t2085;
  double t2086;
  double t2090;
  double t2091;
  double t2092;
  double t2133;
  double t2134;
  double t2135;
  double t2093;
  double t2109;
  double t2113;
  double t2046;
  double t2063;
  double t2064;
  double t2136;
  double t2137;
  double t2146;
  double t2173;
  double t2176;
  double t2185;
  double t2187;
  double t2191;
  double t2192;
  double t2193;
  double t2200;
  double t2201;
  double t2202;
  double t2203;
  double t2204;
  double t2209;
  double t2210;
  double t2219;
  double t2205;
  double t2206;
  double t2207;
  double t2186;
  double t2188;
  double t2189;
  double t2220;
  double t2223;
  double t2224;
  double t1988;
  double t1989;
  double t2008;
  double t2020;
  double t2040;
  double t2129;
  double t2147;
  double t2150;
  double t2083;
  double t2167;
  double t2168;
  double t2169;
  double t2179;
  double t2180;
  double t2181;
  double t2182;
  double t2183;
  double t2208;
  double t2231;
  double t2232;
  double t2194;
  double t2244;
  double t2246;
  double t2248;
  double t2283;
  double t2284;
  double t2285;
  double t2286;
  double t2287;
  double t2288;
  double t2289;
  double t2290;
  double t1901;
  double t1973;
  double t1985;
  double t2044;
  double t2261;
  double t2262;
  double t2158;
  double t2303;
  double t2304;
  double t2305;
  double t2257;
  double t2258;
  double t2259;
  double t2300;
  double t2301;
  double t2302;
  double t2306;
  double t2307;
  double t2253;
  double t2254;
  double t2255;
  double t2256;
  double t2319;
  double t2320;
  double t2292;
  double t2293;
  double t2294;
  double t2295;
  double t2338;
  double t2339;
  double t2340;
  double t2341;
  double t2342;
  double t2343;
  double t2344;
  double t2297;
  double t2298;
  double t2299;
  double t2311;
  double t2313;
  double t2314;
  double t2315;
  double t2355;
  double t2356;
  double t2357;
  double t2321;
  double t2323;
  double t2324;
  double t2325;
  double t2326;
  double t2327;
  double t2328;
  double t2375;
  double t2376;
  double t2377;
  double t2378;
  double t2379;
  double t2380;
  double t2381;
  double t2382;
  double t2175;
  double t2177;
  double t2178;
  double t2184;
  double t2276;
  double t2277;
  double t2240;
  double t2395;
  double t2396;
  double t2397;
  double t2270;
  double t2273;
  double t2274;
  double t2392;
  double t2393;
  double t2394;
  double t2398;
  double t2399;
  double t2264;
  double t2265;
  double t2266;
  double t2267;
  double t2411;
  double t2412;
  double t2384;
  double t2385;
  double t2386;
  double t2387;
  double t2430;
  double t2431;
  double t2432;
  double t2433;
  double t2434;
  double t2435;
  double t2436;
  double t2389;
  double t2390;
  double t2391;
  double t2403;
  double t2405;
  double t2406;
  double t2407;
  double t2447;
  double t2448;
  double t2449;
  double t2413;
  double t2415;
  double t2416;
  double t2417;
  double t2418;
  double t2419;
  double t2420;
  t1838 = Cos(var1[2]);
  t1897 = Cos(var1[3]);
  t1956 = Sin(var1[3]);
  t1917 = Sin(var1[2]);
  t2045 = Cos(var1[4]);
  t2062 = Sin(var1[4]);
  t2066 = t1897*t2045;
  t2070 = -1.*t1956*t2062;
  t2082 = t2066 + t2070;
  t2085 = -1.*t2045;
  t2086 = 1. + t2085;
  t2090 = 0.4*t2086;
  t2091 = 0.64*t2045;
  t2092 = t2090 + t2091;
  t2133 = -1.*t2045*t1956;
  t2134 = -1.*t1897*t2062;
  t2135 = t2133 + t2134;
  t2093 = t2092*t1956;
  t2109 = 0.24*t1897*t2062;
  t2113 = t2093 + t2109;
  t2046 = t2045*t1956;
  t2063 = t1897*t2062;
  t2064 = t2046 + t2063;
  t2136 = t1897*t2092;
  t2137 = -0.24*t1956*t2062;
  t2146 = t2136 + t2137;
  t2173 = Cos(var1[5]);
  t2176 = Sin(var1[5]);
  t2185 = Cos(var1[6]);
  t2187 = Sin(var1[6]);
  t2191 = t2173*t2185;
  t2192 = -1.*t2176*t2187;
  t2193 = t2191 + t2192;
  t2200 = -1.*t2185;
  t2201 = 1. + t2200;
  t2202 = 0.4*t2201;
  t2203 = 0.64*t2185;
  t2204 = t2202 + t2203;
  t2209 = -1.*t2185*t2176;
  t2210 = -1.*t2173*t2187;
  t2219 = t2209 + t2210;
  t2205 = t2204*t2176;
  t2206 = 0.24*t2173*t2187;
  t2207 = t2205 + t2206;
  t2186 = t2185*t2176;
  t2188 = t2173*t2187;
  t2189 = t2186 + t2188;
  t2220 = t2173*t2204;
  t2223 = -0.24*t2176*t2187;
  t2224 = t2220 + t2223;
  t1988 = Power(t1897,2);
  t1989 = 0.11*t1988;
  t2008 = Power(t1956,2);
  t2020 = 0.11*t2008;
  t2040 = t1989 + t2020;
  t2129 = -1.*t2113*t2082;
  t2147 = -1.*t2135*t2146;
  t2150 = t2129 + t2147;
  t2083 = -1.*t1917*t2082;
  t2167 = t2113*t2064;
  t2168 = t2082*t2146;
  t2169 = t2167 + t2168;
  t2179 = Power(t2173,2);
  t2180 = 0.11*t2179;
  t2181 = Power(t2176,2);
  t2182 = 0.11*t2181;
  t2183 = t2180 + t2182;
  t2208 = -1.*t2207*t2193;
  t2231 = -1.*t2219*t2224;
  t2232 = t2208 + t2231;
  t2194 = -1.*t1917*t2193;
  t2244 = t2207*t2189;
  t2246 = t2193*t2224;
  t2248 = t2244 + t2246;
  t2283 = -1.*t2092*t1956;
  t2284 = -0.24*t1897*t2062;
  t2285 = t2283 + t2284;
  t2286 = t2285*t2082;
  t2287 = t2113*t2082;
  t2288 = t2135*t2146;
  t2289 = t2064*t2146;
  t2290 = t2286 + t2287 + t2288 + t2289;
  t1901 = -1.*t1838*t1897;
  t1973 = t1917*t1956;
  t1985 = t1901 + t1973;
  t2044 = -1.*t1985*t2040;
  t2261 = t1838*t2135;
  t2262 = t2261 + t2083;
  t2158 = -1.*t1917*t2135;
  t2303 = -1.*t1897*t2045;
  t2304 = t1956*t2062;
  t2305 = t2303 + t2304;
  t2257 = -1.*t1917*t2064;
  t2258 = t1838*t2082;
  t2259 = t2257 + t2258;
  t2300 = -1.*t2135*t2285;
  t2301 = -1.*t2135*t2113;
  t2302 = -1.*t2082*t2146;
  t2306 = -1.*t2146*t2305;
  t2307 = t2300 + t2301 + t2302 + t2306;
  t2253 = -1.*t1897*t1917;
  t2254 = -1.*t1838*t1956;
  t2255 = t2253 + t2254;
  t2256 = -1.*t2255*t2040;
  t2319 = t1917*t2135;
  t2320 = t2319 + t2258;
  t2292 = t2092*t2045;
  t2293 = Power(t2062,2);
  t2294 = 0.24*t2293;
  t2295 = t2292 + t2294;
  t2338 = -0.24*t2045*t1956;
  t2339 = t2338 + t2284;
  t2340 = t2339*t2082;
  t2341 = 0.24*t1897*t2045;
  t2342 = t2341 + t2137;
  t2343 = t2064*t2342;
  t2344 = t2340 + t2287 + t2288 + t2343;
  t2297 = t2092*t2062;
  t2298 = -0.24*t2045*t2062;
  t2299 = t2297 + t2298;
  t2311 = -1.*t2262*t2150;
  t2313 = t1838*t2305;
  t2314 = t2158 + t2313;
  t2315 = -1.*t2169*t2314;
  t2355 = -1.*t2135*t2339;
  t2356 = -1.*t2082*t2342;
  t2357 = t2355 + t2301 + t2356 + t2306;
  t2321 = -1.*t2320*t2150;
  t2323 = t1917*t2305;
  t2324 = t2261 + t2323;
  t2325 = -1.*t2169*t2324;
  t2326 = t1838*t2064;
  t2327 = t1917*t2082;
  t2328 = t2326 + t2327;
  t2375 = -1.*t2204*t2176;
  t2376 = -0.24*t2173*t2187;
  t2377 = t2375 + t2376;
  t2378 = t2377*t2193;
  t2379 = t2207*t2193;
  t2380 = t2219*t2224;
  t2381 = t2189*t2224;
  t2382 = t2378 + t2379 + t2380 + t2381;
  t2175 = -1.*t1838*t2173;
  t2177 = t1917*t2176;
  t2178 = t2175 + t2177;
  t2184 = -1.*t2178*t2183;
  t2276 = t1838*t2219;
  t2277 = t2276 + t2194;
  t2240 = -1.*t1917*t2219;
  t2395 = -1.*t2173*t2185;
  t2396 = t2176*t2187;
  t2397 = t2395 + t2396;
  t2270 = -1.*t1917*t2189;
  t2273 = t1838*t2193;
  t2274 = t2270 + t2273;
  t2392 = -1.*t2219*t2377;
  t2393 = -1.*t2219*t2207;
  t2394 = -1.*t2193*t2224;
  t2398 = -1.*t2224*t2397;
  t2399 = t2392 + t2393 + t2394 + t2398;
  t2264 = -1.*t2173*t1917;
  t2265 = -1.*t1838*t2176;
  t2266 = t2264 + t2265;
  t2267 = -1.*t2266*t2183;
  t2411 = t1917*t2219;
  t2412 = t2411 + t2273;
  t2384 = t2204*t2185;
  t2385 = Power(t2187,2);
  t2386 = 0.24*t2385;
  t2387 = t2384 + t2386;
  t2430 = -0.24*t2185*t2176;
  t2431 = t2430 + t2376;
  t2432 = t2431*t2193;
  t2433 = 0.24*t2173*t2185;
  t2434 = t2433 + t2223;
  t2435 = t2189*t2434;
  t2436 = t2432 + t2379 + t2380 + t2435;
  t2389 = t2204*t2187;
  t2390 = -0.24*t2185*t2187;
  t2391 = t2389 + t2390;
  t2403 = -1.*t2277*t2232;
  t2405 = t1838*t2397;
  t2406 = t2240 + t2405;
  t2407 = -1.*t2248*t2406;
  t2447 = -1.*t2219*t2431;
  t2448 = -1.*t2193*t2434;
  t2449 = t2447 + t2393 + t2448 + t2398;
  t2413 = -1.*t2412*t2232;
  t2415 = t1917*t2397;
  t2416 = t2276 + t2415;
  t2417 = -1.*t2248*t2416;
  t2418 = t1838*t2189;
  t2419 = t1917*t2193;
  t2420 = t2418 + t2419;
  p_output1[0]=0;
  p_output1[1]=0;
  p_output1[2]=(-0.5*(6.72*t1917 + t2256 - 1.*t2150*t2259 - 1.*t2169*t2262 + t2267 - 1.*t2232*t2274 - 1.*t2248*t2277)*var2[0] - 0.5*(6.72*t1838 + t2044 - 1.*(-1.*t1838*t2064 + t2083)*t2150 - 1.*(-1.*t1838*t2082 + t2158)*t2169 + t2184 - 1.*(-1.*t1838*t2189 + t2194)*t2232 - 1.*(-1.*t1838*t2193 + t2240)*t2248)*var2[1])*var2[2];
  p_output1[3]=var2[2]*(-0.5*(t2256 - 1.*t2290*t2320 + t2321 + t2325 - 1.*t2307*t2328)*var2[0] - 0.5*(t2044 - 1.*t2262*t2290 - 1.*t2259*t2307 + t2311 + t2315)*var2[1] - 0.5*(-2.*t2169*t2290 - 2.*t2150*t2307)*var2[2] - 0.5*(-1.*t2290*t2295 - 1.*t2299*t2307)*var2[3] + 0.12*t2290*var2[4]);
  p_output1[4]=var2[2]*(-0.5*(t2321 + t2325 - 1.*t2320*t2344 - 1.*t2328*t2357)*var2[0] - 0.5*(t2311 + t2315 - 1.*t2262*t2344 - 1.*t2259*t2357)*var2[1] - 0.5*(-2.*t2169*t2344 - 2.*t2150*t2357)*var2[2] - 0.5*(-1.*(0.24*t2045*t2062 - 1.*t2062*t2092)*t2169 - 1.*t2150*(-0.24*Power(t2045,2) + t2292) - 1.*t2295*t2344 - 1.*t2299*t2357)*var2[3] + 0.12*t2344*var2[4]);
  p_output1[5]=var2[2]*(-0.5*(t2267 - 1.*t2382*t2412 + t2413 + t2417 - 1.*t2399*t2420)*var2[0] - 0.5*(t2184 - 1.*t2277*t2382 - 1.*t2274*t2399 + t2403 + t2407)*var2[1] - 0.5*(-2.*t2248*t2382 - 2.*t2232*t2399)*var2[2] - 0.5*(-1.*t2382*t2387 - 1.*t2391*t2399)*var2[5] + 0.12*t2382*var2[6]);
  p_output1[6]=var2[2]*(-0.5*(t2413 + t2417 - 1.*t2412*t2436 - 1.*t2420*t2449)*var2[0] - 0.5*(t2403 + t2407 - 1.*t2277*t2436 - 1.*t2274*t2449)*var2[1] - 0.5*(-2.*t2248*t2436 - 2.*t2232*t2449)*var2[2] - 0.5*(-1.*(0.24*t2185*t2187 - 1.*t2187*t2204)*t2248 - 1.*t2232*(-0.24*Power(t2185,2) + t2384) - 1.*t2387*t2436 - 1.*t2391*t2449)*var2[5] + 0.12*t2436*var2[6]);
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

#include "Ce3_vec3_five_link_walker.hh"

namespace SymFunction
{

void Ce3_vec3_five_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
