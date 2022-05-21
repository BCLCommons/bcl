// this is only included, so that eclipse c++ editor can be used and parses information correctly
// this will not be interpreted, when the program is compiled with an OpenCL compiler
#ifndef __OPENCL_VERSION__
 #include "kernel.h"
#endif

//! @brief gives the weighted results based on distance for the knn prediction
//! @param SORTED_DISTANCES sorted distances
//! @param SORTED_INDECES sorted indeces corresponding to the exp result for that distance
//! @param EXPERIMENTAL_RESULTS known results which is the basis for making the prediction
//! @param OUTPUT the output which is the prediction
//! @param KAPPA the number of indeces, distances, results to consider per query point
//! @param ROWS rows
//! @param COLS cols
//! @param XPAD col padding
//! @param YPAD row padding
__kernel void CalculateWeightedResults
(
  __const __global PRECISION *SORTED_DISTANCES,
  __const __global int       *SORTED_INDECES,
  __const __global PRECISION *EXPERIMENTAL_RESULTS,
          __global PRECISION *OUTPUT,
  __const          uint       KAPPA,
  __const          uint       ROWS,
  __const          uint       COLS,
  __const          uint       XPAD,
  __const          uint       YPAD
)
{
  __const uint xindex = get_global_id( 0);

  if( xindex >= COLS - XPAD) return;
  uint ind;
  PRECISION result = 0, dist, exp, inv_dist, inv_dist_sum = 0;
  for( uint value = 0; value < KAPPA; ++value)
  {
    dist = SORTED_DISTANCES[ mad24( COLS, value, xindex)];
    ind = SORTED_INDECES[ mad24( COLS, value, xindex)];
    exp = EXPERIMENTAL_RESULTS[ ind];
    inv_dist = native_recip( dist);
    inv_dist_sum += inv_dist;
    result += inv_dist * exp;
  }
  result = native_divide( result, inv_dist_sum);
  OUTPUT[ xindex] = result;
}
