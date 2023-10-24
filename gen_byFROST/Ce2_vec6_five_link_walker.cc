/*
 * Automatically Generated from Mathematica.
 * Tue 24 Oct 2023 17:04:47 GMT-04:00
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
  double t1440;
  double t1437;
  double t1438;
  double t1441;
  double t1489;
  double t1392;
  double t1439;
  double t1454;
  double t1487;
  double t1510;
  double t1516;
  double t1517;
  double t1518;
  double t1519;
  double t1488;
  double t1490;
  double t1496;
  double t1497;
  double t1498;
  double t1499;
  double t1501;
  double t1507;
  double t1508;
  double t1537;
  double t1538;
  double t1539;
  double t1523;
  double t1524;
  double t1525;
  double t1526;
  double t1527;
  double t1528;
  double t1541;
  double t1542;
  double t1543;
  double t1545;
  double t1546;
  double t1547;
  double t1548;
  double t1549;
  double t1550;
  double t1565;
  double t1566;
  double t1581;
  double t1582;
  double t1583;
  double t1585;
  double t1586;
  double t1587;
  double t1591;
  double t1592;
  double t1593;
  double t1568;
  double t1569;
  double t1570;
  double t1558;
  double t1559;
  double t1560;
  double t1520;
  double t1521;
  double t1522;
  double t1530;
  double t1531;
  double t1532;
  double t1533;
  double t1540;
  double t1562;
  double t1563;
  double t1555;
  double t1556;
  double t1557;
  double t1561;
  double t1564;
  double t1567;
  double t1571;
  double t1572;
  double t1573;
  double t1575;
  double t1576;
  double t1577;
  double t1578;
  double t1579;
  double t1584;
  double t1588;
  double t1589;
  double t1594;
  double t1595;
  double t1596;
  double t1597;
  double t1598;
  double t1600;
  double t1601;
  double t1602;
  double t1604;
  double t1605;
  double t1606;
  double t1607;
  double t1608;
  double t1626;
  double t1627;
  double t1628;
  double t1629;
  double t1630;
  double t1631;
  double t1580;
  double t1590;
  double t1599;
  double t1603;
  double t1609;
  double t1610;
  double t1615;
  double t1616;
  double t1617;
  double t1618;
  double t1509;
  double t1529;
  double t1534;
  double t1535;
  double t1636;
  double t1637;
  double t1638;
  t1440 = Cos(var1[5]);
  t1437 = Cos(var1[6]);
  t1438 = Sin(var1[5]);
  t1441 = Sin(var1[6]);
  t1489 = Sin(var1[2]);
  t1392 = Cos(var1[2]);
  t1439 = -1.*t1437*t1438;
  t1454 = -1.*t1440*t1441;
  t1487 = t1439 + t1454;
  t1510 = -1.*t1437;
  t1516 = 1. + t1510;
  t1517 = 0.4*t1516;
  t1518 = 0.64*t1437;
  t1519 = t1517 + t1518;
  t1488 = t1392*t1487;
  t1490 = -1.*t1440*t1437;
  t1496 = t1438*t1441;
  t1497 = t1490 + t1496;
  t1498 = t1489*t1497;
  t1499 = t1488 + t1498;
  t1501 = -1.*t1440*t1489;
  t1507 = -1.*t1392*t1438;
  t1508 = t1501 + t1507;
  t1537 = t1392*t1440;
  t1538 = -1.*t1489*t1438;
  t1539 = t1537 + t1538;
  t1523 = t1489*t1487;
  t1524 = t1440*t1437;
  t1525 = -1.*t1438*t1441;
  t1526 = t1524 + t1525;
  t1527 = t1392*t1526;
  t1528 = t1523 + t1527;
  t1541 = t1440*t1489;
  t1542 = t1392*t1438;
  t1543 = t1541 + t1542;
  t1545 = t1437*t1438;
  t1546 = t1440*t1441;
  t1547 = t1545 + t1546;
  t1548 = t1392*t1547;
  t1549 = t1489*t1526;
  t1550 = t1548 + t1549;
  t1565 = -1.*t1489*t1526;
  t1566 = t1488 + t1565;
  t1581 = t1519*t1438;
  t1582 = 0.24*t1440*t1441;
  t1583 = t1581 + t1582;
  t1585 = t1440*t1519;
  t1586 = -0.24*t1438*t1441;
  t1587 = t1585 + t1586;
  t1591 = -1.*t1519*t1438;
  t1592 = -0.24*t1440*t1441;
  t1593 = t1591 + t1592;
  t1568 = -1.*t1489*t1487;
  t1569 = t1392*t1497;
  t1570 = t1568 + t1569;
  t1558 = -1.*t1392*t1440;
  t1559 = t1489*t1438;
  t1560 = t1558 + t1559;
  t1520 = t1519*t1441;
  t1521 = -0.24*t1437*t1441;
  t1522 = t1520 + t1521;
  t1530 = t1519*t1437;
  t1531 = Power(t1441,2);
  t1532 = 0.24*t1531;
  t1533 = t1530 + t1532;
  t1540 = 2.*t1508*t1539;
  t1562 = -1.*t1489*t1547;
  t1563 = t1562 + t1527;
  t1555 = Power(t1508,2);
  t1556 = t1508*t1543;
  t1557 = Power(t1539,2);
  t1561 = t1539*t1560;
  t1564 = t1528*t1563;
  t1567 = t1566*t1550;
  t1571 = t1528*t1570;
  t1572 = t1566*t1499;
  t1573 = t1555 + t1556 + t1557 + t1561 + t1564 + t1567 + t1571 + t1572;
  t1575 = Power(t1440,2);
  t1576 = 0.11*t1575;
  t1577 = Power(t1438,2);
  t1578 = 0.11*t1577;
  t1579 = t1576 + t1578;
  t1584 = -1.*t1583*t1526;
  t1588 = -1.*t1487*t1587;
  t1589 = t1584 + t1588;
  t1594 = t1593*t1526;
  t1595 = t1583*t1526;
  t1596 = t1487*t1587;
  t1597 = t1547*t1587;
  t1598 = t1594 + t1595 + t1596 + t1597;
  t1600 = t1583*t1547;
  t1601 = t1526*t1587;
  t1602 = t1600 + t1601;
  t1604 = -1.*t1487*t1593;
  t1605 = -1.*t1487*t1583;
  t1606 = -1.*t1526*t1587;
  t1607 = -1.*t1587*t1497;
  t1608 = t1604 + t1605 + t1606 + t1607;
  t1626 = t1560*t1579;
  t1627 = t1566*t1589;
  t1628 = t1566*t1598;
  t1629 = t1602*t1570;
  t1630 = t1563*t1608;
  t1631 = t1626 + t1627 + t1628 + t1629 + t1630;
  t1580 = t1508*t1579;
  t1590 = t1528*t1589;
  t1599 = t1528*t1598;
  t1603 = t1602*t1499;
  t1609 = t1550*t1608;
  t1610 = t1580 + t1590 + t1599 + t1603 + t1609;
  t1615 = 0.11*t1560;
  t1616 = t1522*t1566;
  t1617 = t1533*t1570;
  t1618 = t1615 + t1616 + t1617;
  t1509 = 0.11*t1508;
  t1529 = t1522*t1528;
  t1534 = t1533*t1499;
  t1535 = t1509 + t1529 + t1534;
  t1636 = t1533*t1598;
  t1637 = t1522*t1608;
  t1638 = t1636 + t1637;
  p_output1[0]=var2[5]*(-0.5*(2.*t1499*t1528 + t1540 + 2.*t1539*t1543 + 2.*t1528*t1550)*var2[0] - 0.5*t1573*var2[1] - 0.5*t1610*var2[2] - 0.5*t1535*var2[5] - 0.12*t1499*var2[6]);
  p_output1[1]=var2[5]*(-0.5*t1573*var2[0] - 0.5*(t1540 + 2.*t1508*t1560 + 2.*t1563*t1566 + 2.*t1566*t1570)*var2[1] - 0.5*t1631*var2[2] - 0.5*t1618*var2[5] - 0.12*t1570*var2[6]);
  p_output1[2]=var2[5]*(-0.5*t1610*var2[0] - 0.5*t1631*var2[1] - 0.5*(2.*t1598*t1602 + 2.*t1589*t1608)*var2[2] - 0.5*t1638*var2[5] - 0.12*t1598*var2[6]);
  p_output1[3]=0;
  p_output1[4]=0;
  p_output1[5]=(-0.5*t1535*var2[0] - 0.5*t1618*var2[1] - 0.5*t1638*var2[2])*var2[5];
  p_output1[6]=(-0.12*t1499*var2[0] - 0.12*t1570*var2[1] - 0.12*t1598*var2[2])*var2[5];
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

#include "Ce2_vec6_five_link_walker.hh"

namespace SymFunction
{

void Ce2_vec6_five_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
