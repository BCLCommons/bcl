// this is only included, so that eclipse c++ editor can be used and parses information correctly
// this will not be interpreted, when the program is compiled with an OpenCL compiler
#ifndef __OPENCL_VERSION__
 #include "kernel.h"
#endif

//! @brief calculates sigmoid of all elements in matrix except for those added as padding
//! @param  OUTPUT sigmoid output matrix
//! @param  INPUT  input matrix
//! @param  ROWS   rows in input matrix
//! @param  COLS   cols in input matrix
//! @param  XPAD   padding added to cols
//! @param  YPAD   padding added to rows
__kernel void SigmoidKernel
(
          __global PRECISION *OUTPUT,
  __const __global PRECISION *INPUT,
  __const          uint        ROWS,
  __const          uint        COLS,
  __const          uint        XPAD,
  __const          uint        YPAD
)
{
  const uint block_size = get_local_size( 0);
  uint xindex = get_global_id( 0);
  uint yindex = get_global_id( 1);
  uint size_x = COLS - XPAD;
  uint size_y = ROWS - YPAD;
  uint index = yindex * COLS + xindex;
  if( xindex < size_x && yindex < size_y)
  {
    PRECISION temp;
    temp = INPUT[ index];
    OUTPUT[ index] = native_recip( 1 + native_exp( -temp));
  }
}

//! @brief calculates the error in the hidden-to-output
//! @param  OUTPUT output of error calculation
//! @param  INPUT_T target values matrix
//! @param  INPUT_Z predicted values
//! @param  ROWS   rows in input matrix
//! @param  COLS   cols in input matrix
//! @param  XPAD   padding added to cols
//! @param  YPAD   padding added to rows
__kernel void ErrorKKernel
(
          __global PRECISION *OUTPUT,
  __const __global PRECISION *INPUT_T,
  __const __global PRECISION *INPUT_Z,
  __const uint ROWS,
  __const uint COLS,
  __const uint XPAD,
  __const uint YPAD
)
{
  const uint block_size = get_local_size( 0);
  uint xindex = get_global_id( 0);
  uint yindex = get_global_id( 1);
  uint size_x = COLS - XPAD;
  uint size_y = ROWS - YPAD;
  uint index = yindex * COLS + xindex;
  if( xindex < size_x && yindex < size_y)
  {
    PRECISION T, Z;
    T = INPUT_T[ index];
    Z = INPUT_Z[ index];
    OUTPUT[ index] = ( T - Z) * (( Z * ( 1 - Z)));
  }
}

//! @brief  calculates the input-to-hidden errors
//! @param  OUTPUT calculated errors output
//! @param  INPUT_I propagated errors from hidden-to-output errors
//! @param  INPUT_Y output from hidden layer
//! @param  ROWS   rows in input matrix
//! @param  COLS   cols in input matrix
//! @param  XPAD   padding added to cols
//! @param  YPAD   padding added to rows
__kernel void ErrorJKernel
(
          __global PRECISION *OUTPUT,
  __const __global PRECISION *INPUT_I,
  __const __global PRECISION *INPUT_Y,
  __const          uint        ROWS,
  __const          uint        COLS,
  __const          uint        XPAD,
  __const          uint        YPAD
)
{
  const uint block_size = get_local_size( 0);
  uint xindex = get_global_id( 0);
  uint yindex = get_global_id( 1);
  uint size_x = COLS - XPAD;
  uint size_y = ROWS - YPAD;
  uint index = yindex * COLS + xindex;
  if( xindex < size_x && yindex < size_y)
  {
    PRECISION I, Y;
    I = INPUT_I[ index];
    Y = INPUT_Y[ index];
    OUTPUT[ index] = ( I) * (( Y * ( 1 - Y)));
  }
}

//! @brief calculates the change terms
//! @param  DELTA_OUTPUT the change terms
//! @param  SLOPES_INPUT weight changes slopes
//! @param  ETA the learning rate
//! @param  ALPHA the learning momentum
//! @param  ROWS   rows in input matrix
//! @param  COLS   cols in input matrix
//! @param  XPAD   padding added to cols
//! @param  YPAD   padding added to rows
__kernel void DeltaKernel
(
          __global PRECISION *DELTA_OUTPUT,
  __const __global PRECISION *SLOPES_INPUT,
  __const          PRECISION ETA,
  __const          PRECISION ALPHA,
  __const          uint       ROWS,
  __const          uint       COLS,
  __const          uint       XPAD,
  __const          uint       YPAD
)
{
  uint xindex = get_global_id( 0);
  uint yindex = get_global_id( 1);
  uint size_x = COLS - XPAD;
  uint size_y = ROWS - YPAD;
  uint index = yindex * COLS + xindex;
  if( xindex < size_x && yindex < size_y)
  {
    PRECISION slope, delta, eta, alpha;
    eta = ETA;
    alpha = ALPHA;
    slope = SLOPES_INPUT[ index];
    delta = DELTA_OUTPUT[ index];
    DELTA_OUTPUT[ index] = eta * slope + alpha * delta;
  }
}

//! @brief  elementwise addition kernel
//! @param  INPUT_WEIGHTS the weights which will be added to
//! @param  INPUT_CHANGES the values to be added to the weights
//! @param  ROWS   rows in input matrix
//! @param  COLS   cols in input matrix
//! @param  XPAD   padding added to cols
//! @param  YPAD   padding added to rows
__kernel void AddKernel
(
          __global PRECISION *INPUT_WEIGHTS,
  __const __global PRECISION *INPUT_CHANGES,
  __const          uint        ROWS,
  __const          uint        COLS,
  __const          uint        XPAD,
  __const          uint        YPAD
)
{
  const uint block_size = get_local_size( 0);
  uint xindex = get_global_id( 0);
  uint yindex = get_global_id( 1);
  uint size_x = COLS - XPAD;
  uint size_y = ROWS - YPAD;
  uint index = yindex * COLS + xindex;
  if( xindex < size_x && yindex < size_y)
  {
    PRECISION weights, changes;
    weights = INPUT_WEIGHTS[ index];
    changes = INPUT_CHANGES[ index];
    INPUT_WEIGHTS[ index] = weights + changes;
  }
}

//! @brief  adds bias vector to each row in a matrix
//! @param  INPUT  the input matrix
//! @param  BIAS  the bias vector to be added to the input rows
//! @param  ROWS   rows in input matrix
//! @param  COLS   cols in input matrix
//! @param  XPAD   padding added to cols
//! @param  YPAD   padding added to rows
__kernel void AddBiasKernel
(
          __global PRECISION *INPUT,
  __const __global PRECISION *BIAS,
  __const          uint        ROWS,
  __const          uint        COLS,
  __const          uint        XPAD,
  __const          uint        YPAD
)
{
  const uint block_size = get_local_size( 0);
  uint xindex = get_global_id( 0);
  uint yindex = get_global_id( 1);
  uint x_boundary = COLS - XPAD;
  uint y_boundary = ROWS - YPAD;
  uint index = yindex * COLS + xindex;
  if( xindex < x_boundary && yindex < y_boundary  )
  {
    INPUT[ index] += BIAS[ xindex];
  }
}

//! @brief does a partial column-wise reduction sum
//! @param DATASET input matrix to be reduced
//! @param RESULT_MATRIX the resulting partially reduced output matrix
//! @param NUM_DATA_ITEMS total number of elements in input
//! @param VECTOR_LENGTH number of cols in input
//! @param SHARED the shared memory allocation for the reduction
__kernel void SumBiasColumnsKernel
(
           __global PRECISION* DATASET,
           __global PRECISION* RESULT_MATRIX,
  __const           uint        NUM_DATA_ITEMS,
  __const           uint        VECTOR_LENGTH,
           __local  PRECISION* SHARED
)
{
  uint global_index = get_global_id( 0);
  uint vector_index = get_global_id( 1);
  PRECISION sum = 0;
  // Loop sequentially over chunks of input vector
  while( global_index < NUM_DATA_ITEMS)
  {
    PRECISION element = DATASET[ global_index * VECTOR_LENGTH + vector_index];
    sum += element;
    global_index += get_global_size( 0);
  }
  // Perform parallel reduction
  uint local_index = get_local_id( 0);
  SHARED[ local_index] = sum;
  barrier( CLK_LOCAL_MEM_FENCE);
  for( uint offset = get_local_size( 0) / 2; offset > 0; offset = offset / 2)
  {
    if( local_index < offset)
    {
      SHARED[ local_index] += SHARED[ local_index + offset];
    }
    barrier( CLK_LOCAL_MEM_FENCE);
  }
  if( local_index == 0)
  {
    RESULT_MATRIX[ vector_index * get_num_groups( 0) + get_group_id( 0)] = SHARED[ 0];
  }
}

__kernel void ResilientUpdate // RPROP
(
          __global PRECISION *WEIGHTS,
          __global PRECISION *CHANGES,
          __global PRECISION *SLOPES,
          __global PRECISION *PREVSLOPES,
  __const          uint       ROWS,
  __const          uint       COLS,
  __const          uint       XPAD,
  __const          uint       YPAD,
  __const          PRECISION  MIN
)
{
  __const PRECISION increase_factor = 1.2;
  __const PRECISION decrease_factor = 0.5;
  __const PRECISION delta_min = 0.0001;
  __const PRECISION delta_max = 50.0;

  __const uint xindex = get_global_id( 0);
  __const uint yindex = get_global_id( 1);
  __const uint index = yindex * COLS + xindex;

  if( xindex >= COLS - XPAD || yindex >= ROWS - YPAD)
    return;

  // computes new change
  PRECISION next_step = CHANGES[ index];
  PRECISION slope = SLOPES[ index];
  PRECISION previous_slope = PREVSLOPES[ index];
  const PRECISION same_sign = previous_slope * slope;

  // if( same_sign > PRECISION( 0))
  if( same_sign > MIN)
  {
    // next_step may not be zero because then the training will stop
    next_step = min( max( next_step, delta_min) * increase_factor, delta_max);
    WEIGHTS[ index] += slope < 0.0f ? -next_step : next_step;
  }
  // else if( same_sign < PRECISION( 0))
  else if( same_sign < -MIN)
  {
    // next_step may not be zero because then the training will stop
    next_step = max( next_step, delta_min) * decrease_factor;

    // an alternative way to update the weights that accounts for slope direction
    // for unknown reasons, we are not using this at present
    // WEIGHTS[ index] += slope < PRECISION( 0.0) ? -next_step : next_step;

    //slope = 0.0;
    //WEIGHTS[ index] += next_step;

    slope = 0.0;
  }
  else
  {
    WEIGHTS[ index] += slope < 0.0f ? -next_step : next_step;
    //next_step = 0.0;
  }

  // update global data arrays
  PREVSLOPES[ index] = slope;
  CHANGES[ index] = next_step;
  SLOPES[ index] = slope;
}
