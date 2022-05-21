// this is only included, so that eclipse c++ editor can be used and parses information correctly
// this will not be interpreted, when the program is compiled with an OpenCL compiler
#ifndef __OPENCL_VERSION__
 #include "kernel.h"
#endif

__kernel void FindMaxGradient
(
          __global PRECISION *MAX,
          __global uint      *MAX_GRAD_INDEX,
  __const __global PRECISION *GRADIENT,
  __const __global PRECISION *LABELS,
  __const __global uint      *STATUS,
  __const          uint       UPPER_BOUND,
  __const          uint       LOWER_BOUND,
          __local  PRECISION *SHARED,
          __local  uint      *IND_SHARED,
  __const          uint       NUM_ELEMENTS
)
{
  uint gid = get_global_id( 0);
  PRECISION accumulator = -INFINITY;

  // Loop sequentially over chunks of input vector
  while( gid < NUM_ELEMENTS)
  {
    PRECISION label = LABELS[ gid];
    PRECISION gradient = GRADIENT[ gid];
    uint status = STATUS[ gid];
    if( label == 1 && status != UPPER_BOUND && -gradient >= accumulator)
    {
      accumulator = -gradient;
    }
    else if( status != LOWER_BOUND && gradient >= accumulator)
    {
      accumulator = gradient;
    }
    gid += get_global_size( 0);
  }

  // Perform parallel reduction
  int local_index = get_local_id( 0);
  SHARED[ local_index] = accumulator;
  IND_SHARED[ local_index] = get_global_id( 0);
  if( get_global_id( 0) >= NUM_ELEMENTS)
  {
    SHARED[ local_index] = -INFINITY;
    IND_SHARED[ local_index] = 0;
  }
  barrier( CLK_LOCAL_MEM_FENCE);

  for( int offset = get_local_size( 0) / 2; offset > 0; offset = offset / 2)
  {
    if( local_index < offset)
    {
      if( SHARED[ local_index + offset] > SHARED[ local_index])
      {
        SHARED[ local_index] = SHARED[ local_index + offset];
        IND_SHARED[ local_index] = IND_SHARED[ local_index + offset];
      }
      else if( SHARED[ local_index + offset] == SHARED[ local_index])
      {
        IND_SHARED[ local_index + offset] > IND_SHARED[ local_index] ? IND_SHARED[ local_index] = IND_SHARED[ local_index + offset] : 0;
      }
    }
    barrier( CLK_LOCAL_MEM_FENCE);
  }
  if( local_index == 0)
  {
    MAX[ get_group_id( 0)] = SHARED[ 0];
    MAX_GRAD_INDEX[ get_group_id( 0)] = IND_SHARED[ 0];
  }
}

__kernel void ComputeInputIKernelMatrix
(
  __const          uint       VECTOR_INDEX,
  __const __global PRECISION *TRAINING_FEATURES,
  __const          uint       WIDTH,
  __const          uint       HEIGHT,
          __global PRECISION *OUTPUT,
          __local  PRECISION *LOCAL,
                   PRECISION  GAMMA
)
{
  // Each work-group computes multiple elements of RESULT
  for( size_t y = get_group_id( 0); y < HEIGHT; y += get_num_groups( 0))
  {
    const __global PRECISION* row = TRAINING_FEATURES + y * WIDTH;
    const __global PRECISION* vec = TRAINING_FEATURES + VECTOR_INDEX * WIDTH;
    // Each work-item accumulates as many products as necessary
    // into local variable sum
    PRECISION sum = 0, tmp, tmp_row, tmp_vec;
#pragma unroll
    for( size_t x = get_local_id( 0); x < WIDTH; x += get_local_size( 0))
    {
      tmp_row = row[ x], tmp_vec = vec[ x];
      tmp = tmp_row - tmp_vec;
      sum += tmp * tmp;
    }

    // Each partial dot product is stored in shared memory
    LOCAL[ get_local_id( 0)] = sum;

    // Perform parallel reduction to add each work-item's
    // partial dot product together
#pragma unroll
    for( size_t stride = get_local_size( 0) / 2; stride > 0; stride /= 2)
    {
      // Synchronize to make sure each work-item is done updating
      // shared memory; this is necessary because work-items read
      // results written by other work-items during each parallel
      // reduction step
      barrier( CLK_LOCAL_MEM_FENCE);

      // Only the first work-items in the work-group add elements
      // together
      if( get_local_id( 0) < stride)
      {
        // Add two elements from the LOCAL array
        // and store the result in LOCAL[index]
        LOCAL[ get_local_id( 0)] += LOCAL[ get_local_id( 0) + stride];
      }
    }
    // Write the result of the reduction to global memory
    if( get_local_id( 0) == 0)
    {
      OUTPUT[ y] = native_exp( -GAMMA * LOCAL[ 0]);
    }

    // Synchronize to make sure each warp is done reading
    // LOCAL before it is overwritten in the next step
    barrier( CLK_LOCAL_MEM_FENCE);
  }
}

__kernel void GetInputIKernelMatrixResultingVector
(
  __const          uint       VECTOR_ID,
  __const __global PRECISION *KERNEL_VECTOR,
  __const __global int       *SIGNS,
  __const          uint       PROBLEM_LENGTH,
          __global PRECISION *OUTPUT
)
{
  const uint real_prob_size = PROBLEM_LENGTH / 2;
  const uint index = get_global_id( 0);
  if( index >= PROBLEM_LENGTH) return;
  const int sign = SIGNS[ VECTOR_ID];
  const int current_sign = SIGNS[ index];
  PRECISION kernel_value;
  if( index >= real_prob_size)
  {
    kernel_value = KERNEL_VECTOR[ index - real_prob_size];
    OUTPUT[ index] = sign * current_sign * kernel_value;
  }
  else
  {
    kernel_value = KERNEL_VECTOR[ index];
    OUTPUT[ index] = sign * current_sign * kernel_value;
  }
}

__kernel void GetGminGmax2
(
  __const __global PRECISION *Q_I,
  __const __global PRECISION *Q_D,
  __const __global PRECISION *GRADIENT,
  __const __global PRECISION *LABELS,
  __const __global uint      *STATUS,
  __const          uint       UPPER_BOUND,
  __const          uint       LOWER_BOUND,
  __const          uint       INDEX_I,
  __const          PRECISION  MAX_GRADIENT,
  __const          PRECISION  EPS_A,
  __const          uint       NUM_ELEMENTS,
          __global PRECISION *GMAX2_OUTPUT,
          __global PRECISION *OBJ_DIFF_MIN_OUTPUT,
          __global uint      *GMIN_INDEX_OUTOUT,
          __local  PRECISION *GMAX2_SHARED,
          __local  PRECISION *OBJ_DIFF_MIN_SHARED,
          __local  uint      *GMIN_INDEX_SHARED
)
{
  uint gid = get_global_id( 0);
  int index = 0;
  PRECISION gmax2_accumulator = -INFINITY;
  PRECISION obj_diff_min_accumulator = INFINITY;
  uint      gmin_index_accumulator = INFINITY;
  PRECISION grad_diff = 0, obj_diff = 0, quad_coef = 0;
  PRECISION max_grad = MAX_GRADIENT;
  PRECISION q_i = Q_I[ gid];
  int greater_equal, less_than;

  // Loop sequentially over chunks of input vector
  while( gid < NUM_ELEMENTS)
  {
    PRECISION label = LABELS[ gid];
    PRECISION gradient = GRADIENT[ gid];
    uint status = STATUS[ gid];
    if( isequal( label, ( PRECISION)( 1)))
    {
      if( status != LOWER_BOUND)
      {
        grad_diff = MAX_GRADIENT + gradient;
        if( isgreaterequal( -gradient, gmax2_accumulator))
        {
          gmax2_accumulator = gradient;
        }
        if( isgreater( grad_diff, ( PRECISION)( 0)))
        {
          quad_coef = Q_I[ INDEX_I] + Q_D[ gid] - 2 * LABELS[ INDEX_I] * q_i;
          if( isgreater( quad_coef, ( PRECISION)( 0)))
          {
            obj_diff = -native_divide( grad_diff * grad_diff, quad_coef);
          }
          else
          {
            obj_diff = -native_divide( grad_diff * grad_diff, EPS_A);
          }
          if( islessequal( obj_diff, obj_diff_min_accumulator))
          {
            gmin_index_accumulator = gid;
            obj_diff_min_accumulator = obj_diff;
          }
        }
      }
    }
    else
    {
      if( status != UPPER_BOUND)
      {
        grad_diff = max_grad - gradient;
        if( isgreaterequal( -gradient, gmax2_accumulator))
        {
          gmax2_accumulator = -gradient;
        }
        if( isgreater( grad_diff, ( PRECISION)( 0)))
        {
          quad_coef = Q_I[ INDEX_I] + Q_D[ gid] + 2 * LABELS[ INDEX_I] * q_i;
          if( isgreater( quad_coef, ( PRECISION)( 0)))
          {
            obj_diff = -native_divide( grad_diff * grad_diff, quad_coef);
          }
          else
          {
            obj_diff = -native_divide( grad_diff * grad_diff, EPS_A);
          }
          if( islessequal( obj_diff, obj_diff_min_accumulator))
          {
            gmin_index_accumulator = gid;
            obj_diff_min_accumulator = obj_diff;
          }
        }
      }
    }
    gid += get_global_size( 0);
  }

  // Perform parallel reduction
  int local_index = get_local_id( 0);
  GMAX2_SHARED[ local_index] = gmax2_accumulator;
  OBJ_DIFF_MIN_SHARED[ local_index] = obj_diff_min_accumulator;
  GMIN_INDEX_SHARED[ local_index] = get_global_id( 0);
  barrier( CLK_LOCAL_MEM_FENCE);

  for( int offset = get_local_size( 0) / 2; offset > 0; offset = offset / 2)
  {
    if( local_index < offset)
    {
      greater_equal = isgreaterequal( GMAX2_SHARED[ local_index + offset], GMAX2_SHARED[ local_index]);
      GMAX2_SHARED[ local_index] = greater_equal ? GMAX2_SHARED[ local_index + offset] : GMAX2_SHARED[ local_index];
      less_than = isless( OBJ_DIFF_MIN_SHARED[ local_index + offset], OBJ_DIFF_MIN_SHARED[ local_index]);
      OBJ_DIFF_MIN_SHARED[ local_index] = less_than ? OBJ_DIFF_MIN_SHARED[ local_index + offset] : OBJ_DIFF_MIN_SHARED[ local_index];
      GMIN_INDEX_SHARED[ local_index] = less_than ? GMIN_INDEX_SHARED[ local_index + offset] : GMIN_INDEX_SHARED[ local_index];
    }
    barrier( CLK_LOCAL_MEM_FENCE);
  }
  if( local_index == 0)
  {
    GMAX2_OUTPUT[ get_group_id( 0)]        = GMAX2_SHARED[ 0];
    OBJ_DIFF_MIN_OUTPUT[ get_group_id( 0)] = OBJ_DIFF_MIN_SHARED[ 0];
    GMIN_INDEX_OUTOUT[ get_group_id( 0)]   = GMIN_INDEX_SHARED[ 0];
  }
}

__kernel void CalculateBias
(
  __const __global uint      *STATUS,
  __const __global PRECISION *GRADIENT,
  __const __global PRECISION *LABELS,
  __const          uint       LOWER_BOUND,
  __const          uint       UPPER_BOUND,
  __const          uint       NUM_ELEMENTS,
          __local  PRECISION *UPPER_SHARED,
          __local  PRECISION *LOWER_SHARED,
          __local  uint      *NR_FREE_SHARED,
          __local  PRECISION *SUM_FREE_SHARED,
          __global PRECISION *UPPER_OUTPUT,
          __global PRECISION *LOWER_OUTPUT,
          __global uint      *NR_FREE_OUTPUT,
          __global PRECISION *SUM_FREE_OUTPUT
)
{
  int gid = get_global_id( 0);
  PRECISION upper_bound_accumulator = INFINITY;
  PRECISION lower_bound_accumulator = -INFINITY;
  uint      nr_free_accumulator = 0;
  PRECISION sum_free_accumulator = 0;

  // Loop sequentially over chunks of input vector
  while( gid < NUM_ELEMENTS)
  {
    PRECISION label = LABELS[ gid];
    PRECISION gradient = GRADIENT[ gid];
    uint status = STATUS[ gid];
    __const PRECISION label_gradient = label * gradient;

    if( status == LOWER_BOUND)
    {
      if( isgreater( label, ( PRECISION)( 0)))
      {
        upper_bound_accumulator = min( upper_bound_accumulator, label_gradient);
      }
      else
      {
        lower_bound_accumulator = max( lower_bound_accumulator, label_gradient);
      }
    }
    else if( status == UPPER_BOUND)
    {
      if( isless( label, ( PRECISION)( 0)))
      {
        upper_bound_accumulator = min( upper_bound_accumulator, label_gradient);
      }
      else
      {
        lower_bound_accumulator = max( lower_bound_accumulator, label_gradient);
      }
    }
    else
    {
      nr_free_accumulator += 1;
      sum_free_accumulator += label_gradient;
    }
    gid += get_global_size( 0);
  }

  // Perform parallel reduction
  int local_index = get_local_id( 0);
  NR_FREE_SHARED [ local_index] = nr_free_accumulator;
  SUM_FREE_SHARED[ local_index] = sum_free_accumulator;
  UPPER_SHARED   [ local_index] = upper_bound_accumulator;
  LOWER_SHARED   [ local_index] = lower_bound_accumulator;
  barrier( CLK_LOCAL_MEM_FENCE);

  for( int offset = get_local_size( 0) / 2; offset > 0; offset = offset / 2)
  {
    if( local_index < offset)
    {
      NR_FREE_SHARED [ local_index] = NR_FREE_SHARED[ local_index + offset] + NR_FREE_SHARED[ local_index];
      SUM_FREE_SHARED[ local_index] = SUM_FREE_SHARED[ local_index + offset] + SUM_FREE_SHARED[ local_index];
      UPPER_SHARED   [ local_index] = min( UPPER_SHARED[ local_index + offset], UPPER_SHARED[ local_index]);
      LOWER_SHARED   [ local_index] = max( LOWER_SHARED[ local_index + offset], LOWER_SHARED[ local_index]);
    }
    barrier( CLK_LOCAL_MEM_FENCE);
  }
  if( local_index == 0)
  {
    NR_FREE_OUTPUT[ get_group_id( 0)]  = NR_FREE_SHARED[ 0];
    SUM_FREE_OUTPUT[ get_group_id( 0)] = SUM_FREE_SHARED[ 0];
    UPPER_OUTPUT[ get_group_id( 0)]    = UPPER_SHARED[ 0];
    LOWER_OUTPUT[ get_group_id( 0)]    = LOWER_SHARED[ 0];
  }
}

__kernel void UpdateAllAlphaStatus
(
  __const          PRECISION  COST,
  __const          uint       UPPER_BOUND,
  __const          uint       LOWER_BOUND,
  __const          uint       FREE,
          __global PRECISION *ALPHAS,
          __global uint      *STATUS,
  __const          uint       NUM_ELEMENTS
)
{
  const uint gid = get_global_id( 0);
  const PRECISION cost = COST;
  if( gid < NUM_ELEMENTS)
  {
    PRECISION alpha = ALPHAS[ gid];
    if( isgreaterequal( alpha, cost))
    {
      STATUS[ gid] = UPPER_BOUND;
    }
    else if( islessequal( alpha, ( PRECISION)( 0)))
    {
      STATUS[ gid] = LOWER_BOUND;
    }
    else
    {
      STATUS[ gid] = FREE;
    }
  }
}

void UpdateAlphaStatus
(
  __const          uint       INDEX,
          __global PRECISION *ALPHAS,
          __global uint      *STATUS,
  __const          PRECISION  COST,
  __const          uint       UPPER_BOUND,
  __const          uint       LOWER_BOUND,
  __const          uint       FREE
)
{
  PRECISION alpha = ALPHAS[ INDEX];
  if( isgreaterequal( alpha, COST))
  {
    STATUS[ INDEX] = UPPER_BOUND;
  }
  else if( islessequal( alpha, ( PRECISION)( 0)))
  {
    STATUS[ INDEX] = LOWER_BOUND;
  }
  else
  {
    STATUS[ INDEX] = FREE;
  }
}

__kernel void UpdateGradient
(
  __const          PRECISION  dALPHA_I,
  __const          PRECISION  dALPHA_J,
          __global PRECISION *GRADIENT,
  __const __global PRECISION *Q_I,
  __const __global PRECISION *Q_J,
  __const          uint       NUM_ELEMENTS
)
{
  __const uint gid = get_global_id( 0);

  if( gid < NUM_ELEMENTS)
  {
    PRECISION gradient = GRADIENT[ gid];
    PRECISION q_i = Q_I[ gid];
    PRECISION q_j = Q_J[ gid];
    GRADIENT[ gid] += q_i * dALPHA_I + (q_j * dALPHA_J);
  }
}


__kernel void UpdateGradientBar
(
  __const          PRECISION  EFFECTIVE_COST,
          __global PRECISION *GRADIENT_BAR,
  __const __global PRECISION *Q_I,
  __const          uint       NUM_ELEMENTS
)
{
  __const uint gid = get_global_id( 0);

  if( gid < NUM_ELEMENTS)
  {
    PRECISION gradient_bar = GRADIENT_BAR[ gid];
    PRECISION q_i = Q_I[ gid];
    GRADIENT_BAR[ gid] = q_i * EFFECTIVE_COST + gradient_bar;
  }
}


__kernel void UpdateGradientWithBias
(
  __const __global PRECISION *GRADIENT_BAR,
  __const __global PRECISION *BIAS,
  __const          uint       NUM_ELEMENTS,
          __global PRECISION *GRADIENT
)
{
  __const uint gid = get_global_id( 0);

  if( gid < NUM_ELEMENTS)
  {
    PRECISION grad_bar = GRADIENT_BAR[ gid];
    PRECISION bias     = BIAS[ gid];

    GRADIENT[ gid] = grad_bar + bias;
  }
}

__kernel void AddAlphaKernelToGradient
(
  __const          uint       ACTIVE_SIZE,
          __global PRECISION *GRADIENT,
  __const __global PRECISION *Q_I,
  __const          PRECISION  ALPHA,
  __const          uint       NUM_ELEMENTS
)
{
  __const uint gid = get_global_id( 0);

  if( gid < NUM_ELEMENTS)
  {
    PRECISION gradient = GRADIENT[ gid + ACTIVE_SIZE];
    PRECISION alpha = ALPHA;
    PRECISION q_i = Q_I[ gid];
    GRADIENT[ gid + ACTIVE_SIZE] = alpha * q_i + gradient;
  }
}


__kernel void AssembleFinalAlphaVector
(
  __const          uint       DATA_SIZE,
  __const __global PRECISION *INPUT_ALPHAS,
          __global PRECISION *FINAL_ALPHAS,
          __global uint      *SV_INDECES
)
{
  __const uint gid = get_global_id( 0);

  // Loop sequentially over chunks of input vector
  if( gid < DATA_SIZE)
  {
    PRECISION alpha_begin = INPUT_ALPHAS[ gid];
    PRECISION alpha_mid   = INPUT_ALPHAS[ gid + DATA_SIZE];

    __const PRECISION alpha_value = alpha_begin - alpha_mid;

    if( isnotequal( alpha_value, ( PRECISION)( 0)))
    {
      FINAL_ALPHAS[ gid] = alpha_value;

      SV_INDECES[ gid] = gid;
    }
    else
    {
      FINAL_ALPHAS[ gid] = 0;
      SV_INDECES[ gid] = 1e5;
    }
  }
}


__kernel void InitializeGradient
(
  __const __global PRECISION *ALPHAS,
  __const __global PRECISION *Q_I,
          __global PRECISION *GRADIENT,
  __const          uint       NUM_ELEMENTS
)
{

  __const uint gid = get_global_id( 0);

  if( gid < NUM_ELEMENTS)
  {
    PRECISION gradient = GRADIENT[ gid];
    PRECISION alpha    = ALPHAS[ gid];
    PRECISION q_i      = Q_I[ gid];

    GRADIENT[ gid] = alpha * q_i + gradient;
  }
}


__kernel void InitializeGradientBar
(
  __const          PRECISION  COST,
  __const __global PRECISION *Q_I,
          __global PRECISION *GRADIENT_BAR,
  __const          uint       NUM_ELEMENTS
)
{

  __const uint gid = get_global_id( 0);

  if( gid < NUM_ELEMENTS)
  {
    PRECISION gradient = GRADIENT_BAR[ gid];
    PRECISION cost     = COST;
    PRECISION q_i      = Q_I[ gid];

    GRADIENT_BAR[ gid] = cost * q_i + gradient;
  }
}

