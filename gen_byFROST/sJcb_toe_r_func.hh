/*
 * Automatically Generated from Mathematica.
 * Tue 24 Oct 2023 17:04:26 GMT-04:00
 */

#ifndef SJCB_TOE_R_FUNC_HH
#define SJCB_TOE_R_FUNC_HH

#ifdef MATLAB_MEX_FILE
// No need for external definitions
#else // MATLAB_MEX_FILE


#include "math2mat.hpp"
#include "mdefs.hpp"

namespace SymFunction
{

  void sJcb_toe_r_func_raw(double *p_output1, const double *var1);

  inline void sJcb_toe_r_func(Eigen::MatrixXd &p_output1, const Eigen::VectorXd &var1)
  {
    // Check
    // - Inputs
    assert_size_matrix(var1, 7, 1);

	
    // - Outputs
    assert_size_matrix(p_output1, 6, 7);


    // set zero the matrix
    p_output1.setZero();


    // Call Subroutine with raw data
    sJcb_toe_r_func_raw(p_output1.data(), var1.data());
    }
  
  
}

#endif // MATLAB_MEX_FILE

#endif // SJCB_TOE_R_FUNC_HH
