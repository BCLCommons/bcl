// this is only included, so that eclipse c++ editor can be used and parses information correctly
// this will not be interpreted, when the program is compiled with an OpenCL compiler
#ifndef __OPENCL_VERSION__
 #include "kernel.h"
#endif

PRECISION RBFKernelDevice
(
  __const __global PRECISION* FEATURE,
  __const __global PRECISION* SUPPORT_VECTOR,
  __const PRECISION  LENGTH,
  __const PRECISION  GAMMA
)
{
  // initialize the sum
  PRECISION sum = 0.0;

  // iterate over feature indeces
#pragma unroll
  for( size_t i = 0; i < LENGTH; ++i)
  {
    // take difference between feature and support vector
    PRECISION tmp = FEATURE[ i] - SUPPORT_VECTOR[ i];

    // square the difference and add to sum
    sum += tmp * tmp;
  }

  // return the exponent
  return native_exp( -GAMMA * sum);
}

PRECISION PolynomialKernelDevice
(
  __const __global PRECISION* FEATURE,
  __const __global PRECISION* SUPPORT_VECTOR,
  __const PRECISION  LENGTH,
  __const PRECISION  EXPONENT,
  __const PRECISION  INHOMOGENEOUS
)
{
  // initialize the sum
  PRECISION sum = 0.0;

  // iterate over feature indeces
#pragma unroll
  for( size_t i = 0; i < LENGTH; ++i)
  {
    // take the product of feature and support vector and add to sum
    sum += FEATURE[ i] * SUPPORT_VECTOR[ i];
  }

  // return the squaretoot
  return native_powr( sum + INHOMOGENEOUS, EXPONENT);
}

//#define ActualKernel( feature, support_vector, length) RBFKernelDevice(        feature, support_vector, length, 2.2)
//#define ActualKernel( feature, support_vector, length) PolynomialKernelDevice( feature, support_vector, length,   2, 1)
#define ActualKernel( feature, support_vector, length) PLACE_HOLDER_FOR_ACTUAL_KERNEL_CALL

//! @brief classify a single vector
//! the grid layout should be nr_support_vectors, 1, 1 (where nr support vectors is also number of alphas)
//! any padding in the support vector matrix (additional support vectors) needs to have an alpha = 0
//! @param FEATURE the feature to classify
//! @param SUPPORT_VECTORS the support vectors
//! @param ALPHAS the alphas (Lagrange multipliers) - same length as the work item size in x
//! @param LENGTH length of feature vector
//! @param RESULT the result of the classification, vector of length 1
__kernel void
ClassifyFeatureSupportVectors
(
  __const __global PRECISION* FEATURE,
  __const __global PRECISION* SUPPORT_VECTORS,
  __const __global PRECISION* ALPHAS,
  __const          uint       LENGTH,
  __const          uint       ROWS,
          __global PRECISION* RESULT
)
{
  if( get_global_id( 0) >= ROWS) return;
  // current thread
  __const size_t global_id = get_global_id( 0);
  __const __global PRECISION* current_vector = SUPPORT_VECTORS + LENGTH * global_id;

  RESULT[ global_id] = ALPHAS[ global_id] *
    ActualKernel( FEATURE, current_vector, LENGTH);
}
