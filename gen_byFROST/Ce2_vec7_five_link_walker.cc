/*
 * Automatically Generated from Mathematica.
 * Tue 24 Oct 2023 17:04:48 GMT-04:00
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
  double t1552;
  double t1536;
  double t1544;
  double t1553;
  double t1612;
  double t1551;
  double t1554;
  double t1574;
  double t1500;
  double t1623;
  double t1624;
  double t1625;
  double t1632;
  double t1633;
  double t1640;
  double t1641;
  double t1642;
  double t1643;
  double t1644;
  double t1645;
  double t1651;
  double t1611;
  double t1613;
  double t1614;
  double t1619;
  double t1620;
  double t1621;
  double t1655;
  double t1656;
  double t1657;
  double t1658;
  double t1659;
  double t1660;
  double t1675;
  double t1676;
  double t1685;
  double t1686;
  double t1687;
  double t1689;
  double t1690;
  double t1691;
  double t1695;
  double t1696;
  double t1697;
  double t1701;
  double t1702;
  double t1678;
  double t1679;
  double t1680;
  double t1652;
  double t1653;
  double t1654;
  double t1672;
  double t1673;
  double t1634;
  double t1635;
  double t1639;
  double t1647;
  double t1648;
  double t1649;
  double t1662;
  double t1663;
  double t1664;
  double t1674;
  double t1677;
  double t1681;
  double t1682;
  double t1683;
  double t1688;
  double t1692;
  double t1693;
  double t1698;
  double t1699;
  double t1700;
  double t1703;
  double t1704;
  double t1706;
  double t1707;
  double t1708;
  double t1710;
  double t1711;
  double t1712;
  double t1713;
  double t1714;
  double t1732;
  double t1733;
  double t1734;
  double t1735;
  double t1736;
  double t1694;
  double t1705;
  double t1709;
  double t1715;
  double t1716;
  double t1721;
  double t1722;
  double t1723;
  double t1724;
  double t1725;
  double t1646;
  double t1650;
  double t1661;
  double t1665;
  double t1666;
  double t1741;
  double t1742;
  double t1743;
  double t1744;
  double t1745;
  t1552 = Cos(var1[5]);
  t1536 = Cos(var1[6]);
  t1544 = Sin(var1[5]);
  t1553 = Sin(var1[6]);
  t1612 = Sin(var1[2]);
  t1551 = -1.*t1536*t1544;
  t1554 = -1.*t1552*t1553;
  t1574 = t1551 + t1554;
  t1500 = Cos(var1[2]);
  t1623 = -1.*t1536;
  t1624 = 1. + t1623;
  t1625 = 0.4*t1624;
  t1632 = 0.64*t1536;
  t1633 = t1625 + t1632;
  t1640 = t1612*t1574;
  t1641 = t1552*t1536;
  t1642 = -1.*t1544*t1553;
  t1643 = t1641 + t1642;
  t1644 = t1500*t1643;
  t1645 = t1640 + t1644;
  t1651 = t1633*t1536;
  t1611 = t1500*t1574;
  t1613 = -1.*t1552*t1536;
  t1614 = t1544*t1553;
  t1619 = t1613 + t1614;
  t1620 = t1612*t1619;
  t1621 = t1611 + t1620;
  t1655 = t1536*t1544;
  t1656 = t1552*t1553;
  t1657 = t1655 + t1656;
  t1658 = t1500*t1657;
  t1659 = t1612*t1643;
  t1660 = t1658 + t1659;
  t1675 = -1.*t1612*t1643;
  t1676 = t1611 + t1675;
  t1685 = t1633*t1544;
  t1686 = 0.24*t1552*t1553;
  t1687 = t1685 + t1686;
  t1689 = t1552*t1633;
  t1690 = -0.24*t1544*t1553;
  t1691 = t1689 + t1690;
  t1695 = -0.24*t1536*t1544;
  t1696 = -0.24*t1552*t1553;
  t1697 = t1695 + t1696;
  t1701 = 0.24*t1552*t1536;
  t1702 = t1701 + t1690;
  t1678 = -1.*t1612*t1574;
  t1679 = t1500*t1619;
  t1680 = t1678 + t1679;
  t1652 = Power(t1536,2);
  t1653 = -0.24*t1652;
  t1654 = t1651 + t1653;
  t1672 = -1.*t1612*t1657;
  t1673 = t1672 + t1644;
  t1634 = t1633*t1553;
  t1635 = -0.24*t1536*t1553;
  t1639 = t1634 + t1635;
  t1647 = -1.*t1633*t1553;
  t1648 = 0.24*t1536*t1553;
  t1649 = t1647 + t1648;
  t1662 = Power(t1553,2);
  t1663 = 0.24*t1662;
  t1664 = t1651 + t1663;
  t1674 = t1645*t1673;
  t1677 = t1676*t1660;
  t1681 = t1645*t1680;
  t1682 = t1676*t1621;
  t1683 = t1674 + t1677 + t1681 + t1682;
  t1688 = -1.*t1687*t1643;
  t1692 = -1.*t1574*t1691;
  t1693 = t1688 + t1692;
  t1698 = t1697*t1643;
  t1699 = t1687*t1643;
  t1700 = t1574*t1691;
  t1703 = t1657*t1702;
  t1704 = t1698 + t1699 + t1700 + t1703;
  t1706 = t1687*t1657;
  t1707 = t1643*t1691;
  t1708 = t1706 + t1707;
  t1710 = -1.*t1574*t1697;
  t1711 = -1.*t1574*t1687;
  t1712 = -1.*t1643*t1702;
  t1713 = -1.*t1691*t1619;
  t1714 = t1710 + t1711 + t1712 + t1713;
  t1732 = t1676*t1693;
  t1733 = t1676*t1704;
  t1734 = t1708*t1680;
  t1735 = t1673*t1714;
  t1736 = t1732 + t1733 + t1734 + t1735;
  t1694 = t1645*t1693;
  t1705 = t1645*t1704;
  t1709 = t1708*t1621;
  t1715 = t1660*t1714;
  t1716 = t1694 + t1705 + t1709 + t1715;
  t1721 = t1654*t1673;
  t1722 = t1639*t1676;
  t1723 = t1649*t1676;
  t1724 = t1664*t1680;
  t1725 = t1721 + t1722 + t1723 + t1724;
  t1646 = t1639*t1645;
  t1650 = t1649*t1645;
  t1661 = t1654*t1660;
  t1665 = t1664*t1621;
  t1666 = t1646 + t1650 + t1661 + t1665;
  t1741 = t1654*t1693;
  t1742 = t1649*t1708;
  t1743 = t1664*t1704;
  t1744 = t1639*t1714;
  t1745 = t1741 + t1742 + t1743 + t1744;
  p_output1[0]=var2[6]*(-0.5*(2.*t1621*t1645 + 2.*t1645*t1660)*var2[0] - 0.5*t1683*var2[1] - 0.5*t1716*var2[2] - 0.5*t1666*var2[5] - 0.12*t1621*var2[6]);
  p_output1[1]=var2[6]*(-0.5*t1683*var2[0] - 0.5*(2.*t1673*t1676 + 2.*t1676*t1680)*var2[1] - 0.5*t1736*var2[2] - 0.5*t1725*var2[5] - 0.12*t1680*var2[6]);
  p_output1[2]=var2[6]*(-0.5*t1716*var2[0] - 0.5*t1736*var2[1] - 0.5*(2.*t1704*t1708 + 2.*t1693*t1714)*var2[2] - 0.5*t1745*var2[5] - 0.12*t1704*var2[6]);
  p_output1[3]=0;
  p_output1[4]=0;
  p_output1[5]=var2[6]*(-0.5*t1666*var2[0] - 0.5*t1725*var2[1] - 0.5*t1745*var2[2] - 0.5*(2.*t1639*t1654 + 2.*t1649*t1664)*var2[5] - 0.12*t1649*var2[6]);
  p_output1[6]=(-0.12*t1621*var2[0] - 0.12*t1680*var2[1] - 0.12*t1704*var2[2] - 0.12*t1649*var2[5])*var2[6];
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

#include "Ce2_vec7_five_link_walker.hh"

namespace SymFunction
{

void Ce2_vec7_five_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
