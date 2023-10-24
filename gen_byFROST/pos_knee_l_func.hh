/*
 * Automatically Generated from Mathematica.
 * Tue 24 Oct 2023 17:04:20 GMT-04:00
 */

#ifndef POS_KNEE_L_FUNC_HH
#define POS_KNEE_L_FUNC_HH

#ifdef MATLAB_MEX_FILE
// No need for external definitions
#else // MATLAB_MEX_FILE


#include "math2mat.hpp"
#include "mdefs.hpp"

namespace SymFunction
{

  void pos_knee_l_func_raw(double *p_output1, const double *var1);

  inline void pos_knee_l_func(Eigen::MatrixXd &p_output1, const Eigen::VectorXd &var1)
  {
    // Check
    // - Inputs
    assert_size_matrix(var1, 7, 1);

	
    // - Outputs
    assert_size_matrix(p_output1, 1, 3);


    // set zero the matrix
    p_output1.setZero();


    // Call Subroutine with raw data
    pos_knee_l_func_raw(p_output1.data(), var1.data());
    }
  
  
}

#endif // MATLAB_MEX_FILE

#endif // POS_KNEE_L_FUNC_HH
