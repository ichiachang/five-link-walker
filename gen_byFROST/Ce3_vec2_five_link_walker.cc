/*
 * Automatically Generated from Mathematica.
 * Tue 24 Oct 2023 17:04:51 GMT-04:00
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
  double t1766;
  double t1759;
  double t1764;
  double t1780;
  double t1806;
  double t1730;
  double t1765;
  double t1799;
  double t1804;
  double t1805;
  double t1817;
  double t1821;
  double t1832;
  double t1836;
  double t1837;
  double t1853;
  double t1854;
  double t1855;
  double t1856;
  double t1857;
  double t1924;
  double t1921;
  double t1922;
  double t1925;
  double t1923;
  double t1934;
  double t1935;
  double t1938;
  double t1939;
  double t1946;
  double t1947;
  double t1952;
  double t1955;
  double t1964;
  double t1965;
  double t1966;
  double t1967;
  double t1968;
  double t1990;
  double t1991;
  double t1992;
  double t1849;
  double t1850;
  double t1851;
  double t1880;
  double t1876;
  double t1877;
  double t1878;
  double t1879;
  double t1881;
  double t2025;
  double t2028;
  double t2029;
  double t1957;
  double t1959;
  double t1961;
  double t1981;
  double t1977;
  double t1978;
  double t1979;
  double t1980;
  double t1982;
  double t1993;
  double t1994;
  double t1995;
  double t2012;
  double t2011;
  double t2019;
  double t1998;
  double t2007;
  double t2030;
  double t2032;
  double t2034;
  double t2042;
  double t2041;
  double t2043;
  double t2038;
  double t2039;
  double t2100;
  double t2101;
  double t2102;
  double t2104;
  double t2105;
  double t2106;
  double t2120;
  double t2121;
  double t2122;
  double t2124;
  double t2125;
  double t2126;
  double t1852;
  double t1873;
  double t1874;
  double t1875;
  double t1858;
  double t1869;
  double t1870;
  double t1871;
  double t2138;
  double t2139;
  double t2140;
  double t2141;
  double t2142;
  double t1996;
  double t1997;
  double t2047;
  double t2048;
  double t2049;
  double t2050;
  double t2051;
  double t2052;
  double t2053;
  double t2054;
  double t2055;
  double t2058;
  double t2061;
  double t2067;
  double t2068;
  double t2069;
  double t2094;
  double t2095;
  double t2096;
  double t2097;
  double t2098;
  double t2099;
  double t2103;
  double t2107;
  double t2108;
  double t2110;
  double t2111;
  double t2112;
  double t2161;
  double t2162;
  double t2163;
  double t2143;
  double t2144;
  double t2145;
  double t2148;
  double t2149;
  double t2152;
  double t2153;
  double t2154;
  double t2155;
  double t2156;
  double t2157;
  double t2160;
  double t2165;
  double t2166;
  double t2170;
  double t2195;
  double t2196;
  double t2172;
  double t2198;
  double t2199;
  double t2174;
  double t1963;
  double t1974;
  double t1975;
  double t1976;
  double t1969;
  double t1970;
  double t1971;
  double t1972;
  double t2211;
  double t2212;
  double t2213;
  double t2214;
  double t2215;
  double t2036;
  double t2037;
  double t2071;
  double t2072;
  double t2073;
  double t2074;
  double t2075;
  double t2076;
  double t2077;
  double t2078;
  double t2079;
  double t2080;
  double t2081;
  double t2087;
  double t2088;
  double t2089;
  double t2114;
  double t2115;
  double t2116;
  double t2117;
  double t2118;
  double t2119;
  double t2123;
  double t2127;
  double t2128;
  double t2130;
  double t2131;
  double t2132;
  double t2234;
  double t2235;
  double t2236;
  double t2216;
  double t2217;
  double t2218;
  double t2221;
  double t2222;
  double t2225;
  double t2226;
  double t2227;
  double t2228;
  double t2229;
  double t2230;
  double t2233;
  double t2238;
  double t2239;
  double t2243;
  double t2268;
  double t2269;
  double t2245;
  double t2271;
  double t2272;
  double t2247;
  t1766 = Cos(var1[3]);
  t1759 = Cos(var1[4]);
  t1764 = Sin(var1[3]);
  t1780 = Sin(var1[4]);
  t1806 = Cos(var1[2]);
  t1730 = Sin(var1[2]);
  t1765 = -1.*t1759*t1764;
  t1799 = -1.*t1766*t1780;
  t1804 = t1765 + t1799;
  t1805 = -1.*t1730*t1804;
  t1817 = t1766*t1759;
  t1821 = -1.*t1764*t1780;
  t1832 = t1817 + t1821;
  t1836 = -1.*t1806*t1832;
  t1837 = t1805 + t1836;
  t1853 = -1.*t1759;
  t1854 = 1. + t1853;
  t1855 = 0.4*t1854;
  t1856 = 0.64*t1759;
  t1857 = t1855 + t1856;
  t1924 = Cos(var1[5]);
  t1921 = Cos(var1[6]);
  t1922 = Sin(var1[5]);
  t1925 = Sin(var1[6]);
  t1923 = -1.*t1921*t1922;
  t1934 = -1.*t1924*t1925;
  t1935 = t1923 + t1934;
  t1938 = -1.*t1730*t1935;
  t1939 = t1924*t1921;
  t1946 = -1.*t1922*t1925;
  t1947 = t1939 + t1946;
  t1952 = -1.*t1806*t1947;
  t1955 = t1938 + t1952;
  t1964 = -1.*t1921;
  t1965 = 1. + t1964;
  t1966 = 0.4*t1965;
  t1967 = 0.64*t1921;
  t1968 = t1966 + t1967;
  t1990 = -1.*t1766*t1730;
  t1991 = -1.*t1806*t1764;
  t1992 = t1990 + t1991;
  t1849 = -1.*t1806*t1766;
  t1850 = t1730*t1764;
  t1851 = t1849 + t1850;
  t1880 = -1.*t1730*t1832;
  t1876 = t1759*t1764;
  t1877 = t1766*t1780;
  t1878 = t1876 + t1877;
  t1879 = -1.*t1806*t1878;
  t1881 = t1879 + t1880;
  t2025 = -1.*t1924*t1730;
  t2028 = -1.*t1806*t1922;
  t2029 = t2025 + t2028;
  t1957 = -1.*t1806*t1924;
  t1959 = t1730*t1922;
  t1961 = t1957 + t1959;
  t1981 = -1.*t1730*t1947;
  t1977 = t1921*t1922;
  t1978 = t1924*t1925;
  t1979 = t1977 + t1978;
  t1980 = -1.*t1806*t1979;
  t1982 = t1980 + t1981;
  t1993 = t1806*t1766;
  t1994 = -1.*t1730*t1764;
  t1995 = t1993 + t1994;
  t2012 = t1806*t1832;
  t2011 = -1.*t1730*t1878;
  t2019 = t2011 + t2012;
  t1998 = t1806*t1804;
  t2007 = t1998 + t1880;
  t2030 = t1806*t1924;
  t2032 = -1.*t1730*t1922;
  t2034 = t2030 + t2032;
  t2042 = t1806*t1947;
  t2041 = -1.*t1730*t1979;
  t2043 = t2041 + t2042;
  t2038 = t1806*t1935;
  t2039 = t2038 + t1981;
  t2100 = t1857*t1764;
  t2101 = 0.24*t1766*t1780;
  t2102 = t2100 + t2101;
  t2104 = t1766*t1857;
  t2105 = -0.24*t1764*t1780;
  t2106 = t2104 + t2105;
  t2120 = t1968*t1922;
  t2121 = 0.24*t1924*t1925;
  t2122 = t2120 + t2121;
  t2124 = t1924*t1968;
  t2125 = -0.24*t1922*t1925;
  t2126 = t2124 + t2125;
  t1852 = -0.11*t1851;
  t1873 = t1857*t1780;
  t1874 = -0.24*t1759*t1780;
  t1875 = t1873 + t1874;
  t1858 = t1857*t1759;
  t1869 = Power(t1780,2);
  t1870 = 0.24*t1869;
  t1871 = t1858 + t1870;
  t2138 = -1.*t1766*t1759;
  t2139 = t1764*t1780;
  t2140 = t2138 + t2139;
  t2141 = t1806*t2140;
  t2142 = t1805 + t2141;
  t1996 = -2.*t1992*t1995;
  t1997 = -2.*t1992*t1851;
  t2047 = Power(t1992,2);
  t2048 = -1.*t2047;
  t2049 = t1766*t1730;
  t2050 = t1806*t1764;
  t2051 = t2049 + t2050;
  t2052 = -1.*t1992*t2051;
  t2053 = Power(t1995,2);
  t2054 = -1.*t2053;
  t2055 = -1.*t1995*t1851;
  t2058 = t1730*t1804;
  t2061 = t2058 + t2012;
  t2067 = t1806*t1878;
  t2068 = t1730*t1832;
  t2069 = t2067 + t2068;
  t2094 = Power(t1766,2);
  t2095 = 0.11*t2094;
  t2096 = Power(t1764,2);
  t2097 = 0.11*t2096;
  t2098 = t2095 + t2097;
  t2099 = -1.*t1851*t2098;
  t2103 = -1.*t2102*t1832;
  t2107 = -1.*t1804*t2106;
  t2108 = t2103 + t2107;
  t2110 = t2102*t1878;
  t2111 = t1832*t2106;
  t2112 = t2110 + t2111;
  t2161 = -1.*t1857*t1764;
  t2162 = -0.24*t1766*t1780;
  t2163 = t2161 + t2162;
  t2143 = 0.12*var2[4]*t2142;
  t2144 = -1.*t1875*t2007;
  t2145 = -1.*t1871*t2142;
  t2148 = -2.*t2019*t2007;
  t2149 = -2.*t2007*t2142;
  t2152 = -1.*t2061*t2019;
  t2153 = -1.*t2007*t2069;
  t2154 = -1.*t2061*t2142;
  t2155 = t1730*t2140;
  t2156 = t1998 + t2155;
  t2157 = -1.*t2007*t2156;
  t2160 = -1.*t2007*t2108;
  t2165 = t2102*t1832;
  t2166 = t1804*t2106;
  t2170 = -1.*t2112*t2142;
  t2195 = -0.24*t1759*t1764;
  t2196 = t2195 + t2162;
  t2172 = -1.*t1804*t2102;
  t2198 = 0.24*t1766*t1759;
  t2199 = t2198 + t2105;
  t2174 = -1.*t2106*t2140;
  t1963 = -0.11*t1961;
  t1974 = t1968*t1925;
  t1975 = -0.24*t1921*t1925;
  t1976 = t1974 + t1975;
  t1969 = t1968*t1921;
  t1970 = Power(t1925,2);
  t1971 = 0.24*t1970;
  t1972 = t1969 + t1971;
  t2211 = -1.*t1924*t1921;
  t2212 = t1922*t1925;
  t2213 = t2211 + t2212;
  t2214 = t1806*t2213;
  t2215 = t1938 + t2214;
  t2036 = -2.*t2029*t2034;
  t2037 = -2.*t2029*t1961;
  t2071 = Power(t2029,2);
  t2072 = -1.*t2071;
  t2073 = t1924*t1730;
  t2074 = t1806*t1922;
  t2075 = t2073 + t2074;
  t2076 = -1.*t2029*t2075;
  t2077 = Power(t2034,2);
  t2078 = -1.*t2077;
  t2079 = -1.*t2034*t1961;
  t2080 = t1730*t1935;
  t2081 = t2080 + t2042;
  t2087 = t1806*t1979;
  t2088 = t1730*t1947;
  t2089 = t2087 + t2088;
  t2114 = Power(t1924,2);
  t2115 = 0.11*t2114;
  t2116 = Power(t1922,2);
  t2117 = 0.11*t2116;
  t2118 = t2115 + t2117;
  t2119 = -1.*t1961*t2118;
  t2123 = -1.*t2122*t1947;
  t2127 = -1.*t1935*t2126;
  t2128 = t2123 + t2127;
  t2130 = t2122*t1979;
  t2131 = t1947*t2126;
  t2132 = t2130 + t2131;
  t2234 = -1.*t1968*t1922;
  t2235 = -0.24*t1924*t1925;
  t2236 = t2234 + t2235;
  t2216 = 0.12*var2[6]*t2215;
  t2217 = -1.*t1976*t2039;
  t2218 = -1.*t1972*t2215;
  t2221 = -2.*t2043*t2039;
  t2222 = -2.*t2039*t2215;
  t2225 = -1.*t2081*t2043;
  t2226 = -1.*t2039*t2089;
  t2227 = -1.*t2081*t2215;
  t2228 = t1730*t2213;
  t2229 = t2038 + t2228;
  t2230 = -1.*t2039*t2229;
  t2233 = -1.*t2039*t2128;
  t2238 = t2122*t1947;
  t2239 = t1935*t2126;
  t2243 = -1.*t2132*t2215;
  t2268 = -0.24*t1921*t1922;
  t2269 = t2268 + t2235;
  t2245 = -1.*t1935*t2122;
  t2271 = 0.24*t1924*t1921;
  t2272 = t2271 + t2125;
  t2247 = -1.*t2126*t2213;
  p_output1[0]=0;
  p_output1[1]=0;
  p_output1[2]=var2[1]*(-0.5*(-1.*Power(t2007,2) - 1.*Power(t2019,2) - 1.*Power(t2039,2) - 1.*Power(t2043,2) + t2048 + t2052 + t2054 + t2055 - 1.*t1837*t2061 - 1.*t1881*t2069 + t2072 + t2076 + t2078 + t2079 - 1.*t1955*t2081 - 1.*t1982*t2089)*var2[0] - 0.5*(t1996 + t1997 - 2.*t1837*t2007 - 2.*t1881*t2019 + t2036 + t2037 - 2.*t1955*t2039 - 2.*t1982*t2043)*var2[1] - 0.5*(6.72*t1806 + t2099 - 1.*t1881*t2108 - 1.*t1837*t2112 + t2119 - 1.*t1982*t2128 - 1.*t1955*t2132)*var2[2] - 0.5*(t1852 - 1.*t1837*t1871 - 1.*t1875*t1881)*var2[3] + 0.12*t1837*var2[4] - 0.5*(t1963 - 1.*t1955*t1972 - 1.*t1976*t1982)*var2[5] + 0.12*t1955*var2[6]);
  p_output1[3]=var2[1]*(t2143 - 0.5*(t2048 + t2052 + t2054 + t2055 + t2152 + t2153 + t2154 + t2157)*var2[0] - 0.5*(t1996 + t1997 + t2148 + t2149)*var2[1] - 0.5*(t2099 + t2160 - 1.*t2007*(t1878*t2106 + t1832*t2163 + t2165 + t2166) + t2170 - 1.*t2019*(-1.*t1832*t2106 - 1.*t1804*t2163 + t2172 + t2174))*var2[2] - 0.5*(t1852 + t2144 + t2145)*var2[3]);
  p_output1[4]=var2[1]*(t2143 - 0.5*(t2152 + t2153 + t2154 + t2157)*var2[0] - 0.5*(t2148 + t2149)*var2[1] - 0.5*(t2160 + t2170 - 1.*t2019*(t2172 + t2174 - 1.*t1804*t2196 - 1.*t1832*t2199) - 1.*t2007*(t2165 + t2166 + t1832*t2196 + t1878*t2199))*var2[2] - 0.5*(-1.*(0.24*t1759*t1780 - 1.*t1780*t1857)*t2007 - 1.*(-0.24*Power(t1759,2) + t1858)*t2019 + t2144 + t2145)*var2[3]);
  p_output1[5]=var2[1]*(t2216 - 0.5*(t2072 + t2076 + t2078 + t2079 + t2225 + t2226 + t2227 + t2230)*var2[0] - 0.5*(t2036 + t2037 + t2221 + t2222)*var2[1] - 0.5*(t2119 + t2233 - 1.*t2039*(t1979*t2126 + t1947*t2236 + t2238 + t2239) + t2243 - 1.*t2043*(-1.*t1947*t2126 - 1.*t1935*t2236 + t2245 + t2247))*var2[2] - 0.5*(t1963 + t2217 + t2218)*var2[5]);
  p_output1[6]=var2[1]*(t2216 - 0.5*(t2225 + t2226 + t2227 + t2230)*var2[0] - 0.5*(t2221 + t2222)*var2[1] - 0.5*(t2233 + t2243 - 1.*t2043*(t2245 + t2247 - 1.*t1935*t2269 - 1.*t1947*t2272) - 1.*t2039*(t2238 + t2239 + t1947*t2269 + t1979*t2272))*var2[2] - 0.5*(-1.*(0.24*t1921*t1925 - 1.*t1925*t1968)*t2039 - 1.*(-0.24*Power(t1921,2) + t1969)*t2043 + t2217 + t2218)*var2[5]);
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

#include "Ce3_vec2_five_link_walker.hh"

namespace SymFunction
{

void Ce3_vec2_five_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
