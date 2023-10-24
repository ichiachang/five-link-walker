/*
 * Automatically Generated from Mathematica.
 * Tue 24 Oct 2023 17:04:45 GMT-04:00
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
  double t1161;
  double t1113;
  double t1127;
  double t1164;
  double t1183;
  double t1072;
  double t1153;
  double t1168;
  double t1176;
  double t1277;
  double t1278;
  double t1279;
  double t1280;
  double t1281;
  double t1179;
  double t1184;
  double t1185;
  double t1219;
  double t1262;
  double t1263;
  double t1265;
  double t1270;
  double t1271;
  double t1307;
  double t1308;
  double t1309;
  double t1285;
  double t1294;
  double t1295;
  double t1296;
  double t1297;
  double t1298;
  double t1311;
  double t1312;
  double t1313;
  double t1315;
  double t1316;
  double t1317;
  double t1318;
  double t1319;
  double t1320;
  double t1335;
  double t1336;
  double t1351;
  double t1352;
  double t1353;
  double t1355;
  double t1356;
  double t1357;
  double t1361;
  double t1362;
  double t1363;
  double t1338;
  double t1339;
  double t1340;
  double t1328;
  double t1329;
  double t1330;
  double t1282;
  double t1283;
  double t1284;
  double t1300;
  double t1301;
  double t1302;
  double t1303;
  double t1310;
  double t1332;
  double t1333;
  double t1325;
  double t1326;
  double t1327;
  double t1331;
  double t1334;
  double t1337;
  double t1341;
  double t1342;
  double t1343;
  double t1345;
  double t1346;
  double t1347;
  double t1348;
  double t1349;
  double t1354;
  double t1358;
  double t1359;
  double t1364;
  double t1365;
  double t1366;
  double t1367;
  double t1368;
  double t1370;
  double t1371;
  double t1372;
  double t1374;
  double t1375;
  double t1376;
  double t1377;
  double t1378;
  double t1396;
  double t1397;
  double t1398;
  double t1399;
  double t1400;
  double t1401;
  double t1350;
  double t1360;
  double t1369;
  double t1373;
  double t1379;
  double t1380;
  double t1385;
  double t1386;
  double t1387;
  double t1388;
  double t1276;
  double t1299;
  double t1304;
  double t1305;
  double t1406;
  double t1407;
  double t1408;
  t1161 = Cos(var1[3]);
  t1113 = Cos(var1[4]);
  t1127 = Sin(var1[3]);
  t1164 = Sin(var1[4]);
  t1183 = Sin(var1[2]);
  t1072 = Cos(var1[2]);
  t1153 = -1.*t1113*t1127;
  t1168 = -1.*t1161*t1164;
  t1176 = t1153 + t1168;
  t1277 = -1.*t1113;
  t1278 = 1. + t1277;
  t1279 = 0.4*t1278;
  t1280 = 0.64*t1113;
  t1281 = t1279 + t1280;
  t1179 = t1072*t1176;
  t1184 = -1.*t1161*t1113;
  t1185 = t1127*t1164;
  t1219 = t1184 + t1185;
  t1262 = t1183*t1219;
  t1263 = t1179 + t1262;
  t1265 = -1.*t1161*t1183;
  t1270 = -1.*t1072*t1127;
  t1271 = t1265 + t1270;
  t1307 = t1072*t1161;
  t1308 = -1.*t1183*t1127;
  t1309 = t1307 + t1308;
  t1285 = t1183*t1176;
  t1294 = t1161*t1113;
  t1295 = -1.*t1127*t1164;
  t1296 = t1294 + t1295;
  t1297 = t1072*t1296;
  t1298 = t1285 + t1297;
  t1311 = t1161*t1183;
  t1312 = t1072*t1127;
  t1313 = t1311 + t1312;
  t1315 = t1113*t1127;
  t1316 = t1161*t1164;
  t1317 = t1315 + t1316;
  t1318 = t1072*t1317;
  t1319 = t1183*t1296;
  t1320 = t1318 + t1319;
  t1335 = -1.*t1183*t1296;
  t1336 = t1179 + t1335;
  t1351 = t1281*t1127;
  t1352 = 0.24*t1161*t1164;
  t1353 = t1351 + t1352;
  t1355 = t1161*t1281;
  t1356 = -0.24*t1127*t1164;
  t1357 = t1355 + t1356;
  t1361 = -1.*t1281*t1127;
  t1362 = -0.24*t1161*t1164;
  t1363 = t1361 + t1362;
  t1338 = -1.*t1183*t1176;
  t1339 = t1072*t1219;
  t1340 = t1338 + t1339;
  t1328 = -1.*t1072*t1161;
  t1329 = t1183*t1127;
  t1330 = t1328 + t1329;
  t1282 = t1281*t1164;
  t1283 = -0.24*t1113*t1164;
  t1284 = t1282 + t1283;
  t1300 = t1281*t1113;
  t1301 = Power(t1164,2);
  t1302 = 0.24*t1301;
  t1303 = t1300 + t1302;
  t1310 = 2.*t1271*t1309;
  t1332 = -1.*t1183*t1317;
  t1333 = t1332 + t1297;
  t1325 = Power(t1271,2);
  t1326 = t1271*t1313;
  t1327 = Power(t1309,2);
  t1331 = t1309*t1330;
  t1334 = t1298*t1333;
  t1337 = t1336*t1320;
  t1341 = t1298*t1340;
  t1342 = t1336*t1263;
  t1343 = t1325 + t1326 + t1327 + t1331 + t1334 + t1337 + t1341 + t1342;
  t1345 = Power(t1161,2);
  t1346 = 0.11*t1345;
  t1347 = Power(t1127,2);
  t1348 = 0.11*t1347;
  t1349 = t1346 + t1348;
  t1354 = -1.*t1353*t1296;
  t1358 = -1.*t1176*t1357;
  t1359 = t1354 + t1358;
  t1364 = t1363*t1296;
  t1365 = t1353*t1296;
  t1366 = t1176*t1357;
  t1367 = t1317*t1357;
  t1368 = t1364 + t1365 + t1366 + t1367;
  t1370 = t1353*t1317;
  t1371 = t1296*t1357;
  t1372 = t1370 + t1371;
  t1374 = -1.*t1176*t1363;
  t1375 = -1.*t1176*t1353;
  t1376 = -1.*t1296*t1357;
  t1377 = -1.*t1357*t1219;
  t1378 = t1374 + t1375 + t1376 + t1377;
  t1396 = t1330*t1349;
  t1397 = t1336*t1359;
  t1398 = t1336*t1368;
  t1399 = t1372*t1340;
  t1400 = t1333*t1378;
  t1401 = t1396 + t1397 + t1398 + t1399 + t1400;
  t1350 = t1271*t1349;
  t1360 = t1298*t1359;
  t1369 = t1298*t1368;
  t1373 = t1372*t1263;
  t1379 = t1320*t1378;
  t1380 = t1350 + t1360 + t1369 + t1373 + t1379;
  t1385 = 0.11*t1330;
  t1386 = t1284*t1336;
  t1387 = t1303*t1340;
  t1388 = t1385 + t1386 + t1387;
  t1276 = 0.11*t1271;
  t1299 = t1284*t1298;
  t1304 = t1303*t1263;
  t1305 = t1276 + t1299 + t1304;
  t1406 = t1303*t1368;
  t1407 = t1284*t1378;
  t1408 = t1406 + t1407;
  p_output1[0]=var2[3]*(-0.5*(2.*t1263*t1298 + t1310 + 2.*t1309*t1313 + 2.*t1298*t1320)*var2[0] - 0.5*t1343*var2[1] - 0.5*t1380*var2[2] - 0.5*t1305*var2[3] - 0.12*t1263*var2[4]);
  p_output1[1]=var2[3]*(-0.5*t1343*var2[0] - 0.5*(t1310 + 2.*t1271*t1330 + 2.*t1333*t1336 + 2.*t1336*t1340)*var2[1] - 0.5*t1401*var2[2] - 0.5*t1388*var2[3] - 0.12*t1340*var2[4]);
  p_output1[2]=var2[3]*(-0.5*t1380*var2[0] - 0.5*t1401*var2[1] - 0.5*(2.*t1368*t1372 + 2.*t1359*t1378)*var2[2] - 0.5*t1408*var2[3] - 0.12*t1368*var2[4]);
  p_output1[3]=(-0.5*t1305*var2[0] - 0.5*t1388*var2[1] - 0.5*t1408*var2[2])*var2[3];
  p_output1[4]=(-0.12*t1263*var2[0] - 0.12*t1340*var2[1] - 0.12*t1368*var2[2])*var2[3];
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

#include "Ce2_vec4_five_link_walker.hh"

namespace SymFunction
{

void Ce2_vec4_five_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
