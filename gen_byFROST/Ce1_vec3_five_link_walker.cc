/*
 * Automatically Generated from Mathematica.
 * Tue 24 Oct 2023 17:04:34 GMT-04:00
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
  double t378;
  double t436;
  double t379;
  double t433;
  double t542;
  double t561;
  double t565;
  double t568;
  double t569;
  double t560;
  double t562;
  double t563;
  double t574;
  double t575;
  double t578;
  double t579;
  double t580;
  double t564;
  double t572;
  double t573;
  double t581;
  double t584;
  double t585;
  double t587;
  double t588;
  double t589;
  double t631;
  double t634;
  double t637;
  double t619;
  double t624;
  double t626;
  double t645;
  double t646;
  double t647;
  double t586;
  double t614;
  double t617;
  double t618;
  double t629;
  double t630;
  double t641;
  double t642;
  double t643;
  double t644;
  double t648;
  double t669;
  double t673;
  double t674;
  double t675;
  double t680;
  double t694;
  double t697;
  double t684;
  double t701;
  double t702;
  double t687;
  double t432;
  double t437;
  double t438;
  double t441;
  double t520;
  double t526;
  double t530;
  double t536;
  double t537;
  double t729;
  double t731;
  double t740;
  double t742;
  double t754;
  double t755;
  double t758;
  double t761;
  double t762;
  double t763;
  double t764;
  double t765;
  double t770;
  double t771;
  double t772;
  double t766;
  double t767;
  double t768;
  double t741;
  double t743;
  double t749;
  double t773;
  double t774;
  double t775;
  double t730;
  double t732;
  double t733;
  double t734;
  double t735;
  double t736;
  double t737;
  double t738;
  double t739;
  double t759;
  double t769;
  double t776;
  double t777;
  double t788;
  double t789;
  double t782;
  double t783;
  double t784;
  double t779;
  double t791;
  double t792;
  double t793;
  double t800;
  double t801;
  double t802;
  double t790;
  double t795;
  double t796;
  double t803;
  double t804;
  double t805;
  double t806;
  double t807;
  double t808;
  double t817;
  double t818;
  double t810;
  double t820;
  double t821;
  double t812;
  double t726;
  double t727;
  double t628;
  double t638;
  double t639;
  double t721;
  double t722;
  double t683;
  double t685;
  double t689;
  double t837;
  double t700;
  double t703;
  double t704;
  double t839;
  double t840;
  double t841;
  double t842;
  double t706;
  double t707;
  double t708;
  double t833;
  double t834;
  double t835;
  double t836;
  double t780;
  double t857;
  double t858;
  double t859;
  double t860;
  double t781;
  double t794;
  double t797;
  double t798;
  double t864;
  double t753;
  double t760;
  double t809;
  double t811;
  double t813;
  double t870;
  double t819;
  double t822;
  double t823;
  double t872;
  double t873;
  double t874;
  double t825;
  double t826;
  double t827;
  double t902;
  double t903;
  double t904;
  double t905;
  double t907;
  double t908;
  double t909;
  double t931;
  double t932;
  double t933;
  double t934;
  double t936;
  double t937;
  double t938;
  t378 = Cos(var1[3]);
  t436 = Sin(var1[3]);
  t379 = Sin(var1[2]);
  t433 = Cos(var1[2]);
  t542 = Cos(var1[4]);
  t561 = Sin(var1[4]);
  t565 = t378*t542;
  t568 = -1.*t436*t561;
  t569 = t565 + t568;
  t560 = -1.*t542*t436;
  t562 = -1.*t378*t561;
  t563 = t560 + t562;
  t574 = -1.*t542;
  t575 = 1. + t574;
  t578 = 0.4*t575;
  t579 = 0.64*t542;
  t580 = t578 + t579;
  t564 = t379*t563;
  t572 = t433*t569;
  t573 = t564 + t572;
  t581 = t580*t436;
  t584 = 0.24*t378*t561;
  t585 = t581 + t584;
  t587 = t378*t580;
  t588 = -0.24*t436*t561;
  t589 = t587 + t588;
  t631 = t542*t436;
  t634 = t378*t561;
  t637 = t631 + t634;
  t619 = -1.*t580*t436;
  t624 = -0.24*t378*t561;
  t626 = t619 + t624;
  t645 = -1.*t378*t542;
  t646 = t436*t561;
  t647 = t645 + t646;
  t586 = -1.*t585*t569;
  t614 = -1.*t563*t589;
  t617 = t586 + t614;
  t618 = t573*t617;
  t629 = t585*t569;
  t630 = t563*t589;
  t641 = t585*t637;
  t642 = t569*t589;
  t643 = t641 + t642;
  t644 = t433*t563;
  t648 = t379*t647;
  t669 = t644 + t648;
  t673 = t643*t669;
  t674 = t433*t637;
  t675 = t379*t569;
  t680 = t674 + t675;
  t694 = -0.24*t542*t436;
  t697 = t694 + t624;
  t684 = -1.*t563*t585;
  t701 = 0.24*t378*t542;
  t702 = t701 + t588;
  t687 = -1.*t589*t647;
  t432 = -1.*t378*t379;
  t437 = -1.*t433*t436;
  t438 = t432 + t437;
  t441 = Power(t378,2);
  t520 = 0.11*t441;
  t526 = Power(t436,2);
  t530 = 0.11*t526;
  t536 = t520 + t530;
  t537 = t438*t536;
  t729 = Cos(var1[5]);
  t731 = Sin(var1[5]);
  t740 = Cos(var1[6]);
  t742 = Sin(var1[6]);
  t754 = t729*t740;
  t755 = -1.*t731*t742;
  t758 = t754 + t755;
  t761 = -1.*t740;
  t762 = 1. + t761;
  t763 = 0.4*t762;
  t764 = 0.64*t740;
  t765 = t763 + t764;
  t770 = -1.*t740*t731;
  t771 = -1.*t729*t742;
  t772 = t770 + t771;
  t766 = t765*t731;
  t767 = 0.24*t729*t742;
  t768 = t766 + t767;
  t741 = t740*t731;
  t743 = t729*t742;
  t749 = t741 + t743;
  t773 = t729*t765;
  t774 = -0.24*t731*t742;
  t775 = t773 + t774;
  t730 = -1.*t729*t379;
  t732 = -1.*t433*t731;
  t733 = t730 + t732;
  t734 = Power(t729,2);
  t735 = 0.11*t734;
  t736 = Power(t731,2);
  t737 = 0.11*t736;
  t738 = t735 + t737;
  t739 = t733*t738;
  t759 = t433*t758;
  t769 = -1.*t768*t758;
  t776 = -1.*t772*t775;
  t777 = t769 + t776;
  t788 = t379*t772;
  t789 = t788 + t759;
  t782 = t768*t749;
  t783 = t758*t775;
  t784 = t782 + t783;
  t779 = t433*t772;
  t791 = -1.*t765*t731;
  t792 = -0.24*t729*t742;
  t793 = t791 + t792;
  t800 = -1.*t729*t740;
  t801 = t731*t742;
  t802 = t800 + t801;
  t790 = t789*t777;
  t795 = t768*t758;
  t796 = t772*t775;
  t803 = t379*t802;
  t804 = t779 + t803;
  t805 = t784*t804;
  t806 = t433*t749;
  t807 = t379*t758;
  t808 = t806 + t807;
  t817 = -0.24*t740*t731;
  t818 = t817 + t792;
  t810 = -1.*t772*t768;
  t820 = 0.24*t729*t740;
  t821 = t820 + t774;
  t812 = -1.*t775*t802;
  t726 = -1.*t379*t569;
  t727 = t644 + t726;
  t628 = t626*t569;
  t638 = t637*t589;
  t639 = t628 + t629 + t630 + t638;
  t721 = -1.*t379*t637;
  t722 = t721 + t572;
  t683 = -1.*t563*t626;
  t685 = -1.*t569*t589;
  t689 = t683 + t684 + t685 + t687;
  t837 = t727*t617;
  t700 = t697*t569;
  t703 = t637*t702;
  t704 = t700 + t629 + t630 + t703;
  t839 = -1.*t379*t563;
  t840 = t433*t647;
  t841 = t839 + t840;
  t842 = t643*t841;
  t706 = -1.*t563*t697;
  t707 = -1.*t569*t702;
  t708 = t706 + t684 + t707 + t687;
  t833 = -1.*t433*t378;
  t834 = t379*t436;
  t835 = t833 + t834;
  t836 = t835*t536;
  t780 = -1.*t379*t758;
  t857 = -1.*t433*t729;
  t858 = t379*t731;
  t859 = t857 + t858;
  t860 = t859*t738;
  t781 = t779 + t780;
  t794 = t793*t758;
  t797 = t749*t775;
  t798 = t794 + t795 + t796 + t797;
  t864 = -1.*t379*t772;
  t753 = -1.*t379*t749;
  t760 = t753 + t759;
  t809 = -1.*t772*t793;
  t811 = -1.*t758*t775;
  t813 = t809 + t810 + t811 + t812;
  t870 = t781*t777;
  t819 = t818*t758;
  t822 = t749*t821;
  t823 = t819 + t795 + t796 + t822;
  t872 = t433*t802;
  t873 = t864 + t872;
  t874 = t784*t873;
  t825 = -1.*t772*t818;
  t826 = -1.*t758*t821;
  t827 = t825 + t810 + t826 + t812;
  t902 = t580*t542;
  t903 = Power(t561,2);
  t904 = 0.24*t903;
  t905 = t902 + t904;
  t907 = t580*t561;
  t908 = -0.24*t542*t561;
  t909 = t907 + t908;
  t931 = t765*t740;
  t932 = Power(t742,2);
  t933 = 0.24*t932;
  t934 = t931 + t933;
  t936 = t765*t742;
  t937 = -0.24*t740*t742;
  t938 = t936 + t937;
  p_output1[0]=var2[2]*(-0.5*(-6.72*t379 + t537 + t617*t722 + t643*t727 + t739 + t760*t777 + t781*t784)*var2[2] - 0.5*(t537 + t618 + t573*t639 + t673 + t680*t689)*var2[3] - 0.5*(t618 + t673 + t573*t704 + t680*t708)*var2[4] - 0.5*(t739 + t790 + t789*t798 + t805 + t808*t813)*var2[5] - 0.5*(t790 + t805 + t789*t823 + t808*t827)*var2[6]);
  p_output1[1]=var2[2]*(-0.5*(-6.72*t433 + t617*(-1.*t433*t637 + t726) + t777*(-1.*t433*t749 + t780) + t836 + t643*(-1.*t433*t569 + t839) + t860 + t784*(-1.*t433*t758 + t864))*var2[2] - 0.5*(t689*t722 + t639*t727 + t836 + t837 + t842)*var2[3] - 0.5*(t708*t722 + t704*t727 + t837 + t842)*var2[4] - 0.5*(t781*t798 + t760*t813 + t860 + t870 + t874)*var2[5] - 0.5*(t781*t823 + t760*t827 + t870 + t874)*var2[6]);
  p_output1[2]=var2[2]*(-0.5*(2.*t639*t643 + 2.*t617*t689)*var2[3] - 0.5*(2.*t643*t704 + 2.*t617*t708)*var2[4] - 0.5*(2.*t784*t798 + 2.*t777*t813)*var2[5] - 0.5*(2.*t784*t823 + 2.*t777*t827)*var2[6]);
  p_output1[3]=var2[2]*(-0.5*(t639*t905 + t689*t909)*var2[3] - 0.5*((0.24*t542*t561 - 1.*t561*t580)*t643 + t617*(-0.24*Power(t542,2) + t902) + t704*t905 + t708*t909)*var2[4]);
  p_output1[4]=var2[2]*(-0.12*t639*var2[3] - 0.12*t704*var2[4]);
  p_output1[5]=var2[2]*(-0.5*(t798*t934 + t813*t938)*var2[5] - 0.5*((0.24*t740*t742 - 1.*t742*t765)*t784 + t777*(-0.24*Power(t740,2) + t931) + t823*t934 + t827*t938)*var2[6]);
  p_output1[6]=var2[2]*(-0.12*t798*var2[5] - 0.12*t823*var2[6]);
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

#include "Ce1_vec3_five_link_walker.hh"

namespace SymFunction
{

void Ce1_vec3_five_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
