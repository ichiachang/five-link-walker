/*
 * Automatically Generated from Mathematica.
 * Tue 24 Oct 2023 17:04:46 GMT-04:00
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
  double t1322;
  double t1306;
  double t1314;
  double t1323;
  double t1382;
  double t1321;
  double t1324;
  double t1344;
  double t1264;
  double t1393;
  double t1394;
  double t1395;
  double t1402;
  double t1403;
  double t1410;
  double t1411;
  double t1412;
  double t1413;
  double t1414;
  double t1415;
  double t1421;
  double t1381;
  double t1383;
  double t1384;
  double t1389;
  double t1390;
  double t1391;
  double t1425;
  double t1426;
  double t1427;
  double t1428;
  double t1429;
  double t1430;
  double t1445;
  double t1446;
  double t1455;
  double t1456;
  double t1457;
  double t1459;
  double t1460;
  double t1461;
  double t1465;
  double t1466;
  double t1467;
  double t1471;
  double t1472;
  double t1448;
  double t1449;
  double t1450;
  double t1422;
  double t1423;
  double t1424;
  double t1442;
  double t1443;
  double t1404;
  double t1405;
  double t1409;
  double t1417;
  double t1418;
  double t1419;
  double t1432;
  double t1433;
  double t1434;
  double t1444;
  double t1447;
  double t1451;
  double t1452;
  double t1453;
  double t1458;
  double t1462;
  double t1463;
  double t1468;
  double t1469;
  double t1470;
  double t1473;
  double t1474;
  double t1476;
  double t1477;
  double t1478;
  double t1480;
  double t1481;
  double t1482;
  double t1483;
  double t1484;
  double t1502;
  double t1503;
  double t1504;
  double t1505;
  double t1506;
  double t1464;
  double t1475;
  double t1479;
  double t1485;
  double t1486;
  double t1491;
  double t1492;
  double t1493;
  double t1494;
  double t1495;
  double t1416;
  double t1420;
  double t1431;
  double t1435;
  double t1436;
  double t1511;
  double t1512;
  double t1513;
  double t1514;
  double t1515;
  t1322 = Cos(var1[3]);
  t1306 = Cos(var1[4]);
  t1314 = Sin(var1[3]);
  t1323 = Sin(var1[4]);
  t1382 = Sin(var1[2]);
  t1321 = -1.*t1306*t1314;
  t1324 = -1.*t1322*t1323;
  t1344 = t1321 + t1324;
  t1264 = Cos(var1[2]);
  t1393 = -1.*t1306;
  t1394 = 1. + t1393;
  t1395 = 0.4*t1394;
  t1402 = 0.64*t1306;
  t1403 = t1395 + t1402;
  t1410 = t1382*t1344;
  t1411 = t1322*t1306;
  t1412 = -1.*t1314*t1323;
  t1413 = t1411 + t1412;
  t1414 = t1264*t1413;
  t1415 = t1410 + t1414;
  t1421 = t1403*t1306;
  t1381 = t1264*t1344;
  t1383 = -1.*t1322*t1306;
  t1384 = t1314*t1323;
  t1389 = t1383 + t1384;
  t1390 = t1382*t1389;
  t1391 = t1381 + t1390;
  t1425 = t1306*t1314;
  t1426 = t1322*t1323;
  t1427 = t1425 + t1426;
  t1428 = t1264*t1427;
  t1429 = t1382*t1413;
  t1430 = t1428 + t1429;
  t1445 = -1.*t1382*t1413;
  t1446 = t1381 + t1445;
  t1455 = t1403*t1314;
  t1456 = 0.24*t1322*t1323;
  t1457 = t1455 + t1456;
  t1459 = t1322*t1403;
  t1460 = -0.24*t1314*t1323;
  t1461 = t1459 + t1460;
  t1465 = -0.24*t1306*t1314;
  t1466 = -0.24*t1322*t1323;
  t1467 = t1465 + t1466;
  t1471 = 0.24*t1322*t1306;
  t1472 = t1471 + t1460;
  t1448 = -1.*t1382*t1344;
  t1449 = t1264*t1389;
  t1450 = t1448 + t1449;
  t1422 = Power(t1306,2);
  t1423 = -0.24*t1422;
  t1424 = t1421 + t1423;
  t1442 = -1.*t1382*t1427;
  t1443 = t1442 + t1414;
  t1404 = t1403*t1323;
  t1405 = -0.24*t1306*t1323;
  t1409 = t1404 + t1405;
  t1417 = -1.*t1403*t1323;
  t1418 = 0.24*t1306*t1323;
  t1419 = t1417 + t1418;
  t1432 = Power(t1323,2);
  t1433 = 0.24*t1432;
  t1434 = t1421 + t1433;
  t1444 = t1415*t1443;
  t1447 = t1446*t1430;
  t1451 = t1415*t1450;
  t1452 = t1446*t1391;
  t1453 = t1444 + t1447 + t1451 + t1452;
  t1458 = -1.*t1457*t1413;
  t1462 = -1.*t1344*t1461;
  t1463 = t1458 + t1462;
  t1468 = t1467*t1413;
  t1469 = t1457*t1413;
  t1470 = t1344*t1461;
  t1473 = t1427*t1472;
  t1474 = t1468 + t1469 + t1470 + t1473;
  t1476 = t1457*t1427;
  t1477 = t1413*t1461;
  t1478 = t1476 + t1477;
  t1480 = -1.*t1344*t1467;
  t1481 = -1.*t1344*t1457;
  t1482 = -1.*t1413*t1472;
  t1483 = -1.*t1461*t1389;
  t1484 = t1480 + t1481 + t1482 + t1483;
  t1502 = t1446*t1463;
  t1503 = t1446*t1474;
  t1504 = t1478*t1450;
  t1505 = t1443*t1484;
  t1506 = t1502 + t1503 + t1504 + t1505;
  t1464 = t1415*t1463;
  t1475 = t1415*t1474;
  t1479 = t1478*t1391;
  t1485 = t1430*t1484;
  t1486 = t1464 + t1475 + t1479 + t1485;
  t1491 = t1424*t1443;
  t1492 = t1409*t1446;
  t1493 = t1419*t1446;
  t1494 = t1434*t1450;
  t1495 = t1491 + t1492 + t1493 + t1494;
  t1416 = t1409*t1415;
  t1420 = t1419*t1415;
  t1431 = t1424*t1430;
  t1435 = t1434*t1391;
  t1436 = t1416 + t1420 + t1431 + t1435;
  t1511 = t1424*t1463;
  t1512 = t1419*t1478;
  t1513 = t1434*t1474;
  t1514 = t1409*t1484;
  t1515 = t1511 + t1512 + t1513 + t1514;
  p_output1[0]=var2[4]*(-0.5*(2.*t1391*t1415 + 2.*t1415*t1430)*var2[0] - 0.5*t1453*var2[1] - 0.5*t1486*var2[2] - 0.5*t1436*var2[3] - 0.12*t1391*var2[4]);
  p_output1[1]=var2[4]*(-0.5*t1453*var2[0] - 0.5*(2.*t1443*t1446 + 2.*t1446*t1450)*var2[1] - 0.5*t1506*var2[2] - 0.5*t1495*var2[3] - 0.12*t1450*var2[4]);
  p_output1[2]=var2[4]*(-0.5*t1486*var2[0] - 0.5*t1506*var2[1] - 0.5*(2.*t1474*t1478 + 2.*t1463*t1484)*var2[2] - 0.5*t1515*var2[3] - 0.12*t1474*var2[4]);
  p_output1[3]=var2[4]*(-0.5*t1436*var2[0] - 0.5*t1495*var2[1] - 0.5*t1515*var2[2] - 0.5*(2.*t1409*t1424 + 2.*t1419*t1434)*var2[3] - 0.12*t1419*var2[4]);
  p_output1[4]=(-0.12*t1391*var2[0] - 0.12*t1450*var2[1] - 0.12*t1474*var2[2] - 0.12*t1419*var2[3])*var2[4];
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

#include "Ce2_vec5_five_link_walker.hh"

namespace SymFunction
{

void Ce2_vec5_five_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
