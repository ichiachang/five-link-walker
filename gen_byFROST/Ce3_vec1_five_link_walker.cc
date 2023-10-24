/*
 * Automatically Generated from Mathematica.
 * Tue 24 Oct 2023 17:04:49 GMT-04:00
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
  double t1670;
  double t1667;
  double t1668;
  double t1671;
  double t1719;
  double t1622;
  double t1720;
  double t1726;
  double t1727;
  double t1740;
  double t1746;
  double t1747;
  double t1748;
  double t1749;
  double t1669;
  double t1684;
  double t1717;
  double t1718;
  double t1728;
  double t1729;
  double t1770;
  double t1767;
  double t1768;
  double t1771;
  double t1775;
  double t1776;
  double t1777;
  double t1785;
  double t1786;
  double t1787;
  double t1788;
  double t1789;
  double t1769;
  double t1772;
  double t1773;
  double t1774;
  double t1778;
  double t1779;
  double t1731;
  double t1737;
  double t1738;
  double t1807;
  double t1808;
  double t1809;
  double t1757;
  double t1753;
  double t1754;
  double t1755;
  double t1756;
  double t1758;
  double t1781;
  double t1782;
  double t1783;
  double t1822;
  double t1823;
  double t1824;
  double t1797;
  double t1793;
  double t1794;
  double t1795;
  double t1796;
  double t1798;
  double t1811;
  double t1812;
  double t1813;
  double t1815;
  double t1816;
  double t1818;
  double t1819;
  double t1820;
  double t1826;
  double t1827;
  double t1828;
  double t1830;
  double t1831;
  double t1833;
  double t1834;
  double t1835;
  double t1888;
  double t1889;
  double t1890;
  double t1892;
  double t1893;
  double t1894;
  double t1908;
  double t1909;
  double t1910;
  double t1912;
  double t1913;
  double t1914;
  double t1739;
  double t1750;
  double t1751;
  double t1752;
  double t1760;
  double t1761;
  double t1762;
  double t1763;
  double t1926;
  double t1927;
  double t1928;
  double t1929;
  double t1930;
  double t1810;
  double t1814;
  double t1839;
  double t1840;
  double t1841;
  double t1842;
  double t1843;
  double t1844;
  double t1845;
  double t1846;
  double t1847;
  double t1848;
  double t1882;
  double t1883;
  double t1884;
  double t1885;
  double t1886;
  double t1887;
  double t1891;
  double t1895;
  double t1896;
  double t1898;
  double t1899;
  double t1900;
  double t1949;
  double t1950;
  double t1951;
  double t1931;
  double t1932;
  double t1933;
  double t1936;
  double t1937;
  double t1940;
  double t1941;
  double t1942;
  double t1943;
  double t1944;
  double t1945;
  double t1948;
  double t1953;
  double t1954;
  double t1958;
  double t1983;
  double t1984;
  double t1960;
  double t1986;
  double t1987;
  double t1962;
  double t1784;
  double t1790;
  double t1791;
  double t1792;
  double t1800;
  double t1801;
  double t1802;
  double t1803;
  double t1999;
  double t2000;
  double t2001;
  double t2002;
  double t2003;
  double t1825;
  double t1829;
  double t1859;
  double t1860;
  double t1861;
  double t1862;
  double t1863;
  double t1864;
  double t1865;
  double t1866;
  double t1867;
  double t1868;
  double t1902;
  double t1903;
  double t1904;
  double t1905;
  double t1906;
  double t1907;
  double t1911;
  double t1915;
  double t1916;
  double t1918;
  double t1919;
  double t1920;
  double t2022;
  double t2023;
  double t2024;
  double t2004;
  double t2005;
  double t2006;
  double t2009;
  double t2010;
  double t2013;
  double t2014;
  double t2015;
  double t2016;
  double t2017;
  double t2018;
  double t2021;
  double t2026;
  double t2027;
  double t2031;
  double t2056;
  double t2057;
  double t2033;
  double t2059;
  double t2060;
  double t2035;
  t1670 = Cos(var1[3]);
  t1667 = Cos(var1[4]);
  t1668 = Sin(var1[3]);
  t1671 = Sin(var1[4]);
  t1719 = Sin(var1[2]);
  t1622 = Cos(var1[2]);
  t1720 = t1670*t1667;
  t1726 = -1.*t1668*t1671;
  t1727 = t1720 + t1726;
  t1740 = -1.*t1667;
  t1746 = 1. + t1740;
  t1747 = 0.4*t1746;
  t1748 = 0.64*t1667;
  t1749 = t1747 + t1748;
  t1669 = -1.*t1667*t1668;
  t1684 = -1.*t1670*t1671;
  t1717 = t1669 + t1684;
  t1718 = t1622*t1717;
  t1728 = -1.*t1719*t1727;
  t1729 = t1718 + t1728;
  t1770 = Cos(var1[5]);
  t1767 = Cos(var1[6]);
  t1768 = Sin(var1[5]);
  t1771 = Sin(var1[6]);
  t1775 = t1770*t1767;
  t1776 = -1.*t1768*t1771;
  t1777 = t1775 + t1776;
  t1785 = -1.*t1767;
  t1786 = 1. + t1785;
  t1787 = 0.4*t1786;
  t1788 = 0.64*t1767;
  t1789 = t1787 + t1788;
  t1769 = -1.*t1767*t1768;
  t1772 = -1.*t1770*t1771;
  t1773 = t1769 + t1772;
  t1774 = t1622*t1773;
  t1778 = -1.*t1719*t1777;
  t1779 = t1774 + t1778;
  t1731 = -1.*t1670*t1719;
  t1737 = -1.*t1622*t1668;
  t1738 = t1731 + t1737;
  t1807 = t1622*t1670;
  t1808 = -1.*t1719*t1668;
  t1809 = t1807 + t1808;
  t1757 = t1622*t1727;
  t1753 = t1667*t1668;
  t1754 = t1670*t1671;
  t1755 = t1753 + t1754;
  t1756 = -1.*t1719*t1755;
  t1758 = t1756 + t1757;
  t1781 = -1.*t1770*t1719;
  t1782 = -1.*t1622*t1768;
  t1783 = t1781 + t1782;
  t1822 = t1622*t1770;
  t1823 = -1.*t1719*t1768;
  t1824 = t1822 + t1823;
  t1797 = t1622*t1777;
  t1793 = t1767*t1768;
  t1794 = t1770*t1771;
  t1795 = t1793 + t1794;
  t1796 = -1.*t1719*t1795;
  t1798 = t1796 + t1797;
  t1811 = t1670*t1719;
  t1812 = t1622*t1668;
  t1813 = t1811 + t1812;
  t1815 = t1719*t1717;
  t1816 = t1815 + t1757;
  t1818 = t1622*t1755;
  t1819 = t1719*t1727;
  t1820 = t1818 + t1819;
  t1826 = t1770*t1719;
  t1827 = t1622*t1768;
  t1828 = t1826 + t1827;
  t1830 = t1719*t1773;
  t1831 = t1830 + t1797;
  t1833 = t1622*t1795;
  t1834 = t1719*t1777;
  t1835 = t1833 + t1834;
  t1888 = t1749*t1668;
  t1889 = 0.24*t1670*t1671;
  t1890 = t1888 + t1889;
  t1892 = t1670*t1749;
  t1893 = -0.24*t1668*t1671;
  t1894 = t1892 + t1893;
  t1908 = t1789*t1768;
  t1909 = 0.24*t1770*t1771;
  t1910 = t1908 + t1909;
  t1912 = t1770*t1789;
  t1913 = -0.24*t1768*t1771;
  t1914 = t1912 + t1913;
  t1739 = -0.11*t1738;
  t1750 = t1749*t1671;
  t1751 = -0.24*t1667*t1671;
  t1752 = t1750 + t1751;
  t1760 = t1749*t1667;
  t1761 = Power(t1671,2);
  t1762 = 0.24*t1761;
  t1763 = t1760 + t1762;
  t1926 = -1.*t1670*t1667;
  t1927 = t1668*t1671;
  t1928 = t1926 + t1927;
  t1929 = t1719*t1928;
  t1930 = t1718 + t1929;
  t1810 = -2.*t1738*t1809;
  t1814 = -2.*t1813*t1809;
  t1839 = Power(t1738,2);
  t1840 = -1.*t1839;
  t1841 = -1.*t1738*t1813;
  t1842 = Power(t1809,2);
  t1843 = -1.*t1842;
  t1844 = -1.*t1622*t1670;
  t1845 = t1719*t1668;
  t1846 = t1844 + t1845;
  t1847 = -1.*t1809*t1846;
  t1848 = -1.*t1719*t1717;
  t1882 = Power(t1670,2);
  t1883 = 0.11*t1882;
  t1884 = Power(t1668,2);
  t1885 = 0.11*t1884;
  t1886 = t1883 + t1885;
  t1887 = -1.*t1738*t1886;
  t1891 = -1.*t1890*t1727;
  t1895 = -1.*t1717*t1894;
  t1896 = t1891 + t1895;
  t1898 = t1890*t1755;
  t1899 = t1727*t1894;
  t1900 = t1898 + t1899;
  t1949 = -1.*t1749*t1668;
  t1950 = -0.24*t1670*t1671;
  t1951 = t1949 + t1950;
  t1931 = 0.12*var2[4]*t1930;
  t1932 = -1.*t1752*t1816;
  t1933 = -1.*t1763*t1930;
  t1936 = -2.*t1816*t1820;
  t1937 = -2.*t1816*t1930;
  t1940 = -1.*t1816*t1758;
  t1941 = -1.*t1729*t1820;
  t1942 = t1622*t1928;
  t1943 = t1848 + t1942;
  t1944 = -1.*t1816*t1943;
  t1945 = -1.*t1729*t1930;
  t1948 = -1.*t1816*t1896;
  t1953 = t1890*t1727;
  t1954 = t1717*t1894;
  t1958 = -1.*t1900*t1930;
  t1983 = -0.24*t1667*t1668;
  t1984 = t1983 + t1950;
  t1960 = -1.*t1717*t1890;
  t1986 = 0.24*t1670*t1667;
  t1987 = t1986 + t1893;
  t1962 = -1.*t1894*t1928;
  t1784 = -0.11*t1783;
  t1790 = t1789*t1771;
  t1791 = -0.24*t1767*t1771;
  t1792 = t1790 + t1791;
  t1800 = t1789*t1767;
  t1801 = Power(t1771,2);
  t1802 = 0.24*t1801;
  t1803 = t1800 + t1802;
  t1999 = -1.*t1770*t1767;
  t2000 = t1768*t1771;
  t2001 = t1999 + t2000;
  t2002 = t1719*t2001;
  t2003 = t1774 + t2002;
  t1825 = -2.*t1783*t1824;
  t1829 = -2.*t1828*t1824;
  t1859 = Power(t1783,2);
  t1860 = -1.*t1859;
  t1861 = -1.*t1783*t1828;
  t1862 = Power(t1824,2);
  t1863 = -1.*t1862;
  t1864 = -1.*t1622*t1770;
  t1865 = t1719*t1768;
  t1866 = t1864 + t1865;
  t1867 = -1.*t1824*t1866;
  t1868 = -1.*t1719*t1773;
  t1902 = Power(t1770,2);
  t1903 = 0.11*t1902;
  t1904 = Power(t1768,2);
  t1905 = 0.11*t1904;
  t1906 = t1903 + t1905;
  t1907 = -1.*t1783*t1906;
  t1911 = -1.*t1910*t1777;
  t1915 = -1.*t1773*t1914;
  t1916 = t1911 + t1915;
  t1918 = t1910*t1795;
  t1919 = t1777*t1914;
  t1920 = t1918 + t1919;
  t2022 = -1.*t1789*t1768;
  t2023 = -0.24*t1770*t1771;
  t2024 = t2022 + t2023;
  t2004 = 0.12*var2[6]*t2003;
  t2005 = -1.*t1792*t1831;
  t2006 = -1.*t1803*t2003;
  t2009 = -2.*t1831*t1835;
  t2010 = -2.*t1831*t2003;
  t2013 = -1.*t1831*t1798;
  t2014 = -1.*t1779*t1835;
  t2015 = t1622*t2001;
  t2016 = t1868 + t2015;
  t2017 = -1.*t1831*t2016;
  t2018 = -1.*t1779*t2003;
  t2021 = -1.*t1831*t1916;
  t2026 = t1910*t1777;
  t2027 = t1773*t1914;
  t2031 = -1.*t1920*t2003;
  t2056 = -0.24*t1767*t1768;
  t2057 = t2056 + t2023;
  t2033 = -1.*t1773*t1910;
  t2059 = 0.24*t1770*t1767;
  t2060 = t2059 + t1913;
  t2035 = -1.*t1914*t2001;
  p_output1[0]=0;
  p_output1[1]=0;
  p_output1[2]=var2[0]*(-0.5*(t1810 + t1814 - 2.*t1729*t1816 - 2.*t1758*t1820 + t1825 + t1829 - 2.*t1779*t1831 - 2.*t1798*t1835)*var2[0] - 0.5*(-1.*Power(t1729,2) - 1.*Power(t1758,2) - 1.*Power(t1779,2) - 1.*Power(t1798,2) - 1.*(t1728 - 1.*t1622*t1755)*t1820 - 1.*(t1778 - 1.*t1622*t1795)*t1835 + t1840 + t1841 + t1843 + t1847 - 1.*t1816*(-1.*t1622*t1727 + t1848) + t1860 + t1861 + t1863 + t1867 - 1.*t1831*(-1.*t1622*t1777 + t1868))*var2[1] - 0.5*(6.72*t1719 + t1887 - 1.*t1758*t1896 - 1.*t1729*t1900 + t1907 - 1.*t1798*t1916 - 1.*t1779*t1920)*var2[2] - 0.5*(t1739 - 1.*t1752*t1758 - 1.*t1729*t1763)*var2[3] + 0.12*t1729*var2[4] - 0.5*(t1784 - 1.*t1792*t1798 - 1.*t1779*t1803)*var2[5] + 0.12*t1779*var2[6]);
  p_output1[3]=var2[0]*(t1931 - 0.5*(t1810 + t1814 + t1936 + t1937)*var2[0] - 0.5*(t1840 + t1841 + t1843 + t1847 + t1940 + t1941 + t1944 + t1945)*var2[1] - 0.5*(t1887 + t1948 - 1.*t1816*(t1755*t1894 + t1727*t1951 + t1953 + t1954) + t1958 - 1.*t1820*(-1.*t1727*t1894 - 1.*t1717*t1951 + t1960 + t1962))*var2[2] - 0.5*(t1739 + t1932 + t1933)*var2[3]);
  p_output1[4]=var2[0]*(t1931 - 0.5*(t1936 + t1937)*var2[0] - 0.5*(t1940 + t1941 + t1944 + t1945)*var2[1] - 0.5*(t1948 + t1958 - 1.*t1820*(t1960 + t1962 - 1.*t1717*t1984 - 1.*t1727*t1987) - 1.*t1816*(t1953 + t1954 + t1727*t1984 + t1755*t1987))*var2[2] - 0.5*(-1.*(0.24*t1667*t1671 - 1.*t1671*t1749)*t1816 - 1.*(-0.24*Power(t1667,2) + t1760)*t1820 + t1932 + t1933)*var2[3]);
  p_output1[5]=var2[0]*(t2004 - 0.5*(t1825 + t1829 + t2009 + t2010)*var2[0] - 0.5*(t1860 + t1861 + t1863 + t1867 + t2013 + t2014 + t2017 + t2018)*var2[1] - 0.5*(t1907 + t2021 - 1.*t1831*(t1795*t1914 + t1777*t2024 + t2026 + t2027) + t2031 - 1.*t1835*(-1.*t1777*t1914 - 1.*t1773*t2024 + t2033 + t2035))*var2[2] - 0.5*(t1784 + t2005 + t2006)*var2[5]);
  p_output1[6]=var2[0]*(t2004 - 0.5*(t2009 + t2010)*var2[0] - 0.5*(t2013 + t2014 + t2017 + t2018)*var2[1] - 0.5*(t2021 + t2031 - 1.*t1835*(t2033 + t2035 - 1.*t1773*t2057 - 1.*t1777*t2060) - 1.*t1831*(t2026 + t2027 + t1777*t2057 + t1795*t2060))*var2[2] - 0.5*(-1.*(0.24*t1767*t1771 - 1.*t1771*t1789)*t1831 - 1.*(-0.24*Power(t1767,2) + t1800)*t1835 + t2005 + t2006)*var2[5]);
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

#include "Ce3_vec1_five_link_walker.hh"

namespace SymFunction
{

void Ce3_vec1_five_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
