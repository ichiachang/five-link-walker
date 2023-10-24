/*
 * Automatically Generated from Mathematica.
 * Tue 24 Oct 2023 17:04:33 GMT-04:00
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
  double t210;
  double t176;
  double t199;
  double t220;
  double t175;
  double t294;
  double t306;
  double t307;
  double t308;
  double t309;
  double t201;
  double t282;
  double t292;
  double t311;
  double t335;
  double t336;
  double t293;
  double t310;
  double t347;
  double t348;
  double t349;
  double t363;
  double t364;
  double t365;
  double t380;
  double t381;
  double t382;
  double t418;
  double t420;
  double t422;
  double t344;
  double t345;
  double t346;
  double t350;
  double t351;
  double t352;
  double t353;
  double t362;
  double t366;
  double t367;
  double t368;
  double t369;
  double t370;
  double t371;
  double t383;
  double t408;
  double t411;
  double t412;
  double t413;
  double t423;
  double t424;
  double t425;
  double t428;
  double t431;
  double t442;
  double t463;
  double t467;
  double t468;
  double t469;
  double t483;
  double t485;
  double t486;
  double t496;
  double t498;
  double t497;
  double t499;
  double t500;
  double t502;
  double t503;
  double t504;
  double t516;
  double t521;
  double t522;
  double t523;
  double t528;
  double t515;
  double t519;
  double t524;
  double t525;
  double t527;
  double t529;
  double t533;
  double t534;
  double t535;
  double t501;
  double t549;
  double t552;
  double t553;
  double t474;
  double t477;
  double t478;
  double t479;
  double t481;
  double t487;
  double t488;
  double t491;
  double t494;
  double t495;
  double t547;
  double t548;
  double t554;
  double t555;
  double t556;
  double t557;
  double t558;
  double t559;
  double t566;
  double t567;
  double t570;
  double t571;
  double t434;
  double t435;
  double t439;
  double t440;
  double t509;
  double t514;
  double t531;
  double t532;
  double t576;
  double t577;
  double t582;
  double t583;
  double t596;
  double t597;
  double t598;
  double t599;
  double t600;
  double t601;
  double t602;
  double t603;
  double t605;
  double t606;
  double t607;
  double t611;
  double t612;
  double t613;
  double t604;
  double t608;
  double t609;
  double t610;
  double t615;
  double t616;
  double t620;
  double t621;
  double t622;
  double t623;
  double t632;
  double t633;
  double t625;
  double t635;
  double t636;
  double t627;
  double t590;
  double t591;
  double t592;
  double t593;
  double t594;
  double t595;
  double t655;
  double t656;
  double t657;
  double t658;
  double t659;
  double t660;
  double t661;
  double t662;
  double t664;
  double t665;
  double t666;
  double t649;
  double t650;
  double t651;
  double t652;
  double t653;
  double t654;
  double t663;
  double t667;
  double t668;
  double t670;
  double t671;
  double t672;
  double t677;
  double t678;
  double t679;
  double t676;
  double t681;
  double t682;
  double t686;
  double t695;
  double t696;
  double t688;
  double t698;
  double t699;
  double t690;
  double t711;
  double t717;
  double t718;
  double t719;
  double t712;
  double t713;
  double t714;
  double t715;
  double t723;
  double t724;
  double t744;
  double t750;
  double t751;
  double t752;
  double t745;
  double t746;
  double t747;
  double t748;
  double t756;
  double t757;
  t210 = Cos(var1[3]);
  t176 = Cos(var1[4]);
  t199 = Sin(var1[3]);
  t220 = Sin(var1[4]);
  t175 = Sin(var1[2]);
  t294 = Cos(var1[2]);
  t306 = t210*t176;
  t307 = -1.*t199*t220;
  t308 = t306 + t307;
  t309 = t294*t308;
  t201 = -1.*t176*t199;
  t282 = -1.*t210*t220;
  t292 = t201 + t282;
  t311 = t176*t199;
  t335 = t210*t220;
  t336 = t311 + t335;
  t293 = t175*t292;
  t310 = t293 + t309;
  t347 = t294*t292;
  t348 = -1.*t175*t308;
  t349 = t347 + t348;
  t363 = -1.*t210*t176;
  t364 = t199*t220;
  t365 = t363 + t364;
  t380 = -1.*t210*t175;
  t381 = -1.*t294*t199;
  t382 = t380 + t381;
  t418 = t294*t210;
  t420 = -1.*t175*t199;
  t422 = t418 + t420;
  t344 = -1.*t175*t336;
  t345 = t344 + t309;
  t346 = t310*t345;
  t350 = t294*t336;
  t351 = t175*t308;
  t352 = t350 + t351;
  t353 = t349*t352;
  t362 = -1.*t175*t292;
  t366 = t294*t365;
  t367 = t362 + t366;
  t368 = t310*t367;
  t369 = t175*t365;
  t370 = t347 + t369;
  t371 = t349*t370;
  t383 = Power(t382,2);
  t408 = t210*t175;
  t411 = t294*t199;
  t412 = t408 + t411;
  t413 = t382*t412;
  t423 = Power(t422,2);
  t424 = -1.*t294*t210;
  t425 = t175*t199;
  t428 = t424 + t425;
  t431 = t422*t428;
  t442 = Cos(var1[5]);
  t463 = -1.*t442*t175;
  t467 = Sin(var1[5]);
  t468 = -1.*t294*t467;
  t469 = t463 + t468;
  t483 = t294*t442;
  t485 = -1.*t175*t467;
  t486 = t483 + t485;
  t496 = Cos(var1[6]);
  t498 = Sin(var1[6]);
  t497 = -1.*t496*t467;
  t499 = -1.*t442*t498;
  t500 = t497 + t499;
  t502 = t442*t496;
  t503 = -1.*t467*t498;
  t504 = t502 + t503;
  t516 = t294*t504;
  t521 = t496*t467;
  t522 = t442*t498;
  t523 = t521 + t522;
  t528 = -1.*t175*t504;
  t515 = t175*t500;
  t519 = t515 + t516;
  t524 = -1.*t175*t523;
  t525 = t524 + t516;
  t527 = t294*t500;
  t529 = t527 + t528;
  t533 = t294*t523;
  t534 = t175*t504;
  t535 = t533 + t534;
  t501 = -1.*t175*t500;
  t549 = -1.*t442*t496;
  t552 = t467*t498;
  t553 = t549 + t552;
  t474 = Power(t469,2);
  t477 = t442*t175;
  t478 = t294*t467;
  t479 = t477 + t478;
  t481 = t469*t479;
  t487 = Power(t486,2);
  t488 = -1.*t294*t442;
  t491 = t175*t467;
  t494 = t488 + t491;
  t495 = t486*t494;
  t547 = t519*t525;
  t548 = t529*t535;
  t554 = t294*t553;
  t555 = t501 + t554;
  t556 = t519*t555;
  t557 = t175*t553;
  t558 = t527 + t557;
  t559 = t529*t558;
  t566 = 2.*t345*t349;
  t567 = 2.*t349*t367;
  t570 = 2.*t382*t422;
  t571 = 2.*t382*t428;
  t434 = -1.*t294*t308;
  t435 = t362 + t434;
  t439 = -1.*t294*t336;
  t440 = t439 + t348;
  t509 = -1.*t294*t504;
  t514 = t501 + t509;
  t531 = -1.*t294*t523;
  t532 = t531 + t528;
  t576 = 2.*t469*t486;
  t577 = 2.*t469*t494;
  t582 = 2.*t525*t529;
  t583 = 2.*t529*t555;
  t596 = -1.*t176;
  t597 = 1. + t596;
  t598 = 0.4*t597;
  t599 = 0.64*t176;
  t600 = t598 + t599;
  t601 = t600*t199;
  t602 = 0.24*t210*t220;
  t603 = t601 + t602;
  t605 = t210*t600;
  t606 = -0.24*t199*t220;
  t607 = t605 + t606;
  t611 = -1.*t600*t199;
  t612 = -0.24*t210*t220;
  t613 = t611 + t612;
  t604 = -1.*t603*t308;
  t608 = -1.*t292*t607;
  t609 = t604 + t608;
  t610 = t349*t609;
  t615 = t603*t308;
  t616 = t292*t607;
  t620 = t603*t336;
  t621 = t308*t607;
  t622 = t620 + t621;
  t623 = t622*t367;
  t632 = -0.24*t176*t199;
  t633 = t632 + t612;
  t625 = -1.*t292*t603;
  t635 = 0.24*t210*t176;
  t636 = t635 + t606;
  t627 = -1.*t607*t365;
  t590 = Power(t210,2);
  t591 = 0.11*t590;
  t592 = Power(t199,2);
  t593 = 0.11*t592;
  t594 = t591 + t593;
  t595 = t428*t594;
  t655 = -1.*t496;
  t656 = 1. + t655;
  t657 = 0.4*t656;
  t658 = 0.64*t496;
  t659 = t657 + t658;
  t660 = t659*t467;
  t661 = 0.24*t442*t498;
  t662 = t660 + t661;
  t664 = t442*t659;
  t665 = -0.24*t467*t498;
  t666 = t664 + t665;
  t649 = Power(t442,2);
  t650 = 0.11*t649;
  t651 = Power(t467,2);
  t652 = 0.11*t651;
  t653 = t650 + t652;
  t654 = t494*t653;
  t663 = -1.*t662*t504;
  t667 = -1.*t500*t666;
  t668 = t663 + t667;
  t670 = t662*t523;
  t671 = t504*t666;
  t672 = t670 + t671;
  t677 = -1.*t659*t467;
  t678 = -0.24*t442*t498;
  t679 = t677 + t678;
  t676 = t529*t668;
  t681 = t662*t504;
  t682 = t500*t666;
  t686 = t672*t555;
  t695 = -0.24*t496*t467;
  t696 = t695 + t678;
  t688 = -1.*t500*t662;
  t698 = 0.24*t442*t496;
  t699 = t698 + t665;
  t690 = -1.*t666*t553;
  t711 = 0.11*t428;
  t717 = t600*t220;
  t718 = -0.24*t176*t220;
  t719 = t717 + t718;
  t712 = t600*t176;
  t713 = Power(t220,2);
  t714 = 0.24*t713;
  t715 = t712 + t714;
  t723 = t719*t349;
  t724 = t715*t367;
  t744 = 0.11*t494;
  t750 = t659*t498;
  t751 = -0.24*t496*t498;
  t752 = t750 + t751;
  t745 = t659*t496;
  t746 = Power(t498,2);
  t747 = 0.24*t746;
  t748 = t745 + t747;
  t756 = t752*t529;
  t757 = t748*t555;
  p_output1[0]=var2[1]*(-0.5*(Power(t345,2) + Power(t349,2) + t383 + t413 + t423 + t431 + t310*t435 + t352*t440 + t474 + t481 + t487 + t495 + t514*t519 + Power(t525,2) + Power(t529,2) + t532*t535)*var2[2] - 0.5*(t346 + t353 + t368 + t371 + t383 + t413 + t423 + t431)*var2[3] - 0.5*(t346 + t353 + t368 + t371)*var2[4] - 0.5*(t474 + t481 + t487 + t495 + t547 + t548 + t556 + t559)*var2[5] - 0.5*(t547 + t548 + t556 + t559)*var2[6]);
  p_output1[1]=var2[1]*(-0.5*(2.*t349*t435 + 2.*t345*t440 + 2.*t514*t529 + 2.*t525*t532 + t570 + t571 + t576 + t577)*var2[2] - 0.5*(t566 + t567 + t570 + t571)*var2[3] - 0.5*(t566 + t567)*var2[4] - 0.5*(t576 + t577 + t582 + t583)*var2[5] - 0.5*(t582 + t583)*var2[6]);
  p_output1[2]=var2[1]*(-0.5*(-6.72*t294 + t595 + t440*t609 + t435*t622 + t654 + t532*t668 + t514*t672)*var2[2] - 0.5*(t595 + t610 + t349*(t336*t607 + t308*t613 + t615 + t616) + t623 + t345*(-1.*t308*t607 - 1.*t292*t613 + t625 + t627))*var2[3] - 0.5*(t610 + t623 + t345*(t625 + t627 - 1.*t292*t633 - 1.*t308*t636) + t349*(t615 + t616 + t308*t633 + t336*t636))*var2[4] - 0.5*(t654 + t676 + t529*(t523*t666 + t504*t679 + t681 + t682) + t686 + t525*(-1.*t504*t666 - 1.*t500*t679 + t688 + t690))*var2[5] - 0.5*(t676 + t686 + t525*(t688 + t690 - 1.*t500*t696 - 1.*t504*t699) + t529*(t681 + t682 + t504*t696 + t523*t699))*var2[6]);
  p_output1[3]=var2[1]*(-0.5*(t711 + t435*t715 + t440*t719)*var2[2] - 0.5*(t711 + t723 + t724)*var2[3] - 0.5*(t349*(0.24*t176*t220 - 1.*t220*t600) + t345*(-0.24*Power(t176,2) + t712) + t723 + t724)*var2[4]);
  p_output1[4]=var2[1]*(-0.12*t435*var2[2] - 0.12*t367*var2[3] - 0.12*t367*var2[4]);
  p_output1[5]=var2[1]*(-0.5*(t744 + t514*t748 + t532*t752)*var2[2] - 0.5*(t744 + t756 + t757)*var2[5] - 0.5*(t529*(0.24*t496*t498 - 1.*t498*t659) + t525*(-0.24*Power(t496,2) + t745) + t756 + t757)*var2[6]);
  p_output1[6]=var2[1]*(-0.12*t514*var2[2] - 0.12*t555*var2[5] - 0.12*t555*var2[6]);
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

#include "Ce1_vec2_five_link_walker.hh"

namespace SymFunction
{

void Ce1_vec2_five_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
