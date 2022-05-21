// this is only included, so that eclipse c++ editor can be used and parses information correctly
// this will not be interpreted, when the program is compiled with an OpenCL compiler
#ifndef __OPENCL_VERSION__
 #include "kernel.h"
#endif

//! @brief calculate correlation between two density maps
//! @param EXP_INPUT array with experimental intensities
//! @param SIM_INPUT array with simulated intensities
//! @param CONTOUR_LEVEL intensity below which the exp/sim values are not considered for the calculation
//! @param NUM_ELEMENTS number of elements in the input arrays
//! @param EXP_SUM_OUTPUT     sum of experimental values
//! @param SIM_SUM_OUTPUT     sum of simulated values
//! @param EXP_SIM_SUM_OUTPUT sum of exp and sim values product
//! @param EXP2_SUM_OUTPUT    sum of exp^2
//! @param SIM2_SUM_OUTPUT    sum of sim^2
//! @param COUNT_OUTPUT       number of elements considered, that were above the contour level
//! @param SH_EXP_SUM         local memory
//! @param SH_SIM_SUM         local memory
//! @param SH_EXP_SIM_SUM     local memory
//! @param SH_EXP2_SUM        local memory
//! @param SH_SIM2_SUM        local memory
//! @param SH_COUNT_SUM       local memory
__kernel void DensityCorrelation
(
  const __global PRECISION* EXP_INPUT,
  const __global PRECISION* SIM_INPUT,
  const          PRECISION  CONTOUR_LEVEL,
  const          uint       NUM_ELEMENTS,
        __global PRECISION* EXP_SUM_OUTPUT,
        __global PRECISION* SIM_SUM_OUTPUT,
        __global PRECISION* EXP_SIM_SUM_OUTPUT,
        __global PRECISION* EXP2_SUM_OUTPUT,
        __global PRECISION* SIM2_SUM_OUTPUT,
        __global int*       COUNT_OUTPUT,
        __local  PRECISION* SH_EXP_SUM,
        __local  PRECISION* SH_SIM_SUM,
        __local  PRECISION* SH_EXP_SIM_SUM,
        __local  PRECISION* SH_EXP2_SUM,
        __local  PRECISION* SH_SIM2_SUM,
        __local  int*       SH_COUNT_SUM
)
{
  // acquire position in grid and sizes
  __const size_t block_size = get_local_size( 0);
  __const size_t thread_id  = get_local_id( 0);
  __const size_t grid_size  = block_size * get_num_groups( 0);

  // start index for that thread
  size_t index            = get_group_id( 0) * block_size + get_local_id( 0);

  // Initialize elements in shared memory to 0
  SH_EXP_SUM    [ thread_id] = 0;
  SH_SIM_SUM    [ thread_id] = 0;
  SH_EXP_SIM_SUM[ thread_id] = 0;
  SH_EXP2_SUM   [ thread_id] = 0;
  SH_SIM2_SUM   [ thread_id] = 0;
  SH_COUNT_SUM  [ thread_id] = 0;

  // perform first level of reduction,
  // reading from global memory, writing to shared memory
  while( index < NUM_ELEMENTS)
  {
    // simulated intensity
    __const PRECISION tmp_sim = SIM_INPUT[ index];
    if( tmp_sim > CONTOUR_LEVEL)
    {
      // experimental density
      __const PRECISION tmp_exp = EXP_INPUT[ index];

      // add products
      SH_EXP_SUM    [ thread_id] += tmp_exp;
      SH_SIM_SUM    [ thread_id] += tmp_sim;
      SH_EXP_SIM_SUM[ thread_id] += tmp_exp * tmp_sim;
      SH_EXP2_SUM   [ thread_id] += tmp_exp * tmp_exp;
      SH_SIM2_SUM   [ thread_id] += tmp_sim * tmp_sim;
      SH_COUNT_SUM  [ thread_id] += 1;
    }

    // go to next block
    index += grid_size;
  }

  // wait for block
  barrier( CLK_LOCAL_MEM_FENCE);

  // do reduction in shared memory for each group
  // reduce within block
  for( size_t offset = block_size / 2; offset > 0; offset = offset / 2)
  {
    if( thread_id < offset)
    {
      SH_EXP_SUM    [ thread_id] += SH_EXP_SUM    [ thread_id + offset];
      SH_SIM_SUM    [ thread_id] += SH_SIM_SUM    [ thread_id + offset];
      SH_EXP_SIM_SUM[ thread_id] += SH_EXP_SIM_SUM[ thread_id + offset];
      SH_EXP2_SUM   [ thread_id] += SH_EXP2_SUM   [ thread_id + offset];
      SH_SIM2_SUM   [ thread_id] += SH_SIM2_SUM   [ thread_id + offset];
      SH_COUNT_SUM  [ thread_id] += SH_COUNT_SUM  [ thread_id + offset];
      barrier( CLK_LOCAL_MEM_FENCE);
    }
  }

  // copy result to output
  if( thread_id == 0)
  {
    EXP_SUM_OUTPUT    [ get_group_id( 0)] = SH_EXP_SUM    [ 0];
    SIM_SUM_OUTPUT    [ get_group_id( 0)] = SH_SIM_SUM    [ 0];
    EXP_SIM_SUM_OUTPUT[ get_group_id( 0)] = SH_EXP_SIM_SUM[ 0];
    EXP2_SUM_OUTPUT   [ get_group_id( 0)] = SH_EXP2_SUM   [ 0];
    SIM2_SUM_OUTPUT   [ get_group_id( 0)] = SH_SIM2_SUM   [ 0];
    COUNT_OUTPUT      [ get_group_id( 0)] = SH_COUNT_SUM  [ 0];
  }
}

//! @brief calculate correlation between two density maps optimized for CPU
//! @param EXP_INPUT array with experimental intensities
//! @param SIM_INPUT array with simulated intensities
//! @param CONTOUR_LEVEL intensity below which the exp/sim values are not considered for the calculation
//! @param NUM_ELEMENTS number of elements in the input arrays
//! @param EXP_SUM_OUTPUT     sum of experimental values
//! @param SIM_SUM_OUTPUT     sum of simulated values
//! @param EXP_SIM_SUM_OUTPUT sum of exp and sim values product
//! @param EXP2_SUM_OUTPUT    sum of exp^2
//! @param SIM2_SUM_OUTPUT    sum of sim^2
//! @param COUNT_OUTPUT       number of elements considered, that were above the contour level
__kernel void DensityCorrelationCPU
(
  const __global PRECISION* EXP_INPUT,
  const __global PRECISION* SIM_INPUT,
  const          PRECISION  CONTOUR_LEVEL,
  const          uint       NUM_ELEMENTS,
  const          uint       BLOCK,
        __global PRECISION* EXP_SUM_OUTPUT,
        __global PRECISION* SIM_SUM_OUTPUT,
        __global PRECISION* EXP_SIM_SUM_OUTPUT,
        __global PRECISION* EXP2_SUM_OUTPUT,
        __global PRECISION* SIM2_SUM_OUTPUT,
        __global int*       COUNT_OUTPUT
)
{
  // acquire position in grid and sizes
  unsigned int global_index = get_global_id( 0) * BLOCK;
  int upper_bound = ( get_global_id( 0) + 1) * BLOCK;
  if( upper_bound > NUM_ELEMENTS)
  {
    upper_bound = NUM_ELEMENTS;
  }

  // Initialize elements in shared memory to 0
  PRECISION tmp_exp_sum     = 0;
  PRECISION tmp_sim_sum     = 0;
  PRECISION tmp_exp_sim_sum = 0;
  PRECISION tmp_exp_sum2    = 0;
  PRECISION tmp_sim_sum2    = 0;
  PRECISION tmp_count       = 0;
  PRECISION tmp_sim;
  PRECISION tmp_exp;

  while( global_index < upper_bound)
  {
    // simulated intensity
    tmp_sim = SIM_INPUT[ global_index];
    if( tmp_sim > CONTOUR_LEVEL)
    {
      // experimental density
      tmp_exp      = EXP_INPUT[ global_index];

      // add products
      tmp_exp_sum     += tmp_exp;
      tmp_sim_sum     += tmp_sim;
      tmp_exp_sim_sum += tmp_exp * tmp_sim;
      tmp_exp_sum2    += tmp_exp * tmp_exp;
      tmp_sim_sum2    += tmp_sim * tmp_sim;
      tmp_count       += 1;
    }

    ++global_index;
  }

  EXP_SUM_OUTPUT    [ get_group_id( 0)] = tmp_exp_sum    ;
  SIM_SUM_OUTPUT    [ get_group_id( 0)] = tmp_sim_sum    ;
  EXP_SIM_SUM_OUTPUT[ get_group_id( 0)] = tmp_exp_sim_sum;
  EXP2_SUM_OUTPUT   [ get_group_id( 0)] = tmp_exp_sum2   ;
  SIM2_SUM_OUTPUT   [ get_group_id( 0)] = tmp_sim_sum2   ;
  COUNT_OUTPUT      [ get_group_id( 0)] = tmp_count      ;
}

//! @brief calculate correlation between two density maps of different sizes
//! @param EXP_INPUT array with experimental intensities
//! @param SIM_INPUT array with simulated intensities
//! @param CONTOUR_LEVEL intensity below which the exp/sim values are not considered for the calculation
//! @param NUM_ELEMENTS number of elements in the input arrays
//! @param EXP_SUM_OUTPUT     sum of experimental values
//! @param SIM_SUM_OUTPUT     sum of simulated values
//! @param EXP_SIM_SUM_OUTPUT sum of exp and sim values product
//! @param EXP2_SUM_OUTPUT    sum of exp^2
//! @param SIM2_SUM_OUTPUT    sum of sim^2
//! @param COUNT_OUTPUT       number of elements considered, that were above the contour level
//! @param SH_EXP_SUM         local memory
//! @param SH_SIM_SUM         local memory
//! @param SH_EXP_SIM_SUM     local memory
//! @param SH_EXP2_SUM        local memory
//! @param SH_SIM2_SUM        local memory
//! @param SH_COUNT_SUM       local memory
__kernel void DensityCorrelationOverlap
(
  const __global PRECISION* EXP_INPUT,
  const __global PRECISION* SIM_INPUT,
  const          PRECISION  CONTOUR_LEVEL,
  const          int        EXP_START_X,
  const          int        EXP_START_Y,
  const          int        EXP_START_Z,
  const          uint       EXP_DIM_X,
  const          uint       EXP_DIM_Y,
  const          uint       EXP_DIM_Z,
  const          int        SIM_START_X,
  const          int        SIM_START_Y,
  const          int        SIM_START_Z,
  const          uint       SIM_DIM_X,
  const          uint       SIM_DIM_Y,
  const          uint       SIM_DIM_Z,
  const          uint       EXTENSION_X,
  const          uint       EXTENSION_Y,
  const          uint       EXTENSION_Z,
        __global PRECISION* EXP_SUM_OUTPUT,
        __global PRECISION* SIM_SUM_OUTPUT,
        __global PRECISION* EXP_SIM_SUM_OUTPUT,
        __global PRECISION* EXP2_SUM_OUTPUT,
        __global PRECISION* SIM2_SUM_OUTPUT,
        __global int*       COUNT_OUTPUT,
        __local  PRECISION* SH_EXP_SUM,
        __local  PRECISION* SH_SIM_SUM,
        __local  PRECISION* SH_EXP_SIM_SUM,
        __local  PRECISION* SH_EXP2_SUM,
        __local  PRECISION* SH_SIM2_SUM,
        __local  int*       SH_COUNT_SUM
)
{
  // acquire position in grid and sizes
  __const size_t block_size_x = get_local_size( 0);
  __const size_t block_size_y = get_local_size( 1);
  __const size_t block_size_z = get_local_size( 2);
  __const size_t block_size   = block_size_x * block_size_y * block_size_z;
  __const size_t thread_id_x  = get_local_id( 0);
  __const size_t thread_id_y  = get_local_id( 1);
  __const size_t thread_id_z  = get_local_id( 2);
  __const size_t thread_id    = thread_id_z * block_size_x * block_size_y + thread_id_y * block_size_x + thread_id_x;
  __const size_t grid_size_x  = block_size_x * get_num_groups( 0);
  __const size_t grid_size_y  = block_size_y * get_num_groups( 1);
  __const size_t grid_size_z  = block_size_z * get_num_groups( 2);

  // Initialize elements in shared memory to 0
  SH_EXP_SUM    [ thread_id] = 0;
  SH_SIM_SUM    [ thread_id] = 0;
  SH_EXP_SIM_SUM[ thread_id] = 0;
  SH_EXP2_SUM   [ thread_id] = 0;
  SH_SIM2_SUM   [ thread_id] = 0;
  SH_COUNT_SUM  [ thread_id] = 0;

  // perform first level of reduction,
  // reading from global memory, writing to shared memory
  // start index for that thread
  size_t index_x            = get_group_id( 0) * block_size_x + thread_id_x;
  while( index_x < EXTENSION_X)
  {
    size_t index_y            = get_group_id( 1) * block_size_y + thread_id_y;
    while( index_y < EXTENSION_Y)
    {
      size_t index_z            = get_group_id( 2) * block_size_z + thread_id_z;
      while( index_z < EXTENSION_Z)
      {
        __const size_t index_sim = ( index_z + SIM_START_Z) * SIM_DIM_X * SIM_DIM_Y + ( index_y + SIM_START_Y) * SIM_DIM_X + SIM_START_X + index_x;
        // simulated intensity
        __const PRECISION tmp_sim = SIM_INPUT[ index_sim];
        if( tmp_sim > CONTOUR_LEVEL)
        {
          __const size_t index_exp = ( index_z + EXP_START_Z) * EXP_DIM_X * EXP_DIM_Y + ( index_y + EXP_START_Y) * EXP_DIM_X + EXP_START_X + index_x;
          // experimental density
          __const PRECISION tmp_exp = EXP_INPUT[ index_exp];

          // add products
          SH_EXP_SUM    [ thread_id] += tmp_exp;
          SH_SIM_SUM    [ thread_id] += tmp_sim;
          SH_EXP_SIM_SUM[ thread_id] += tmp_exp * tmp_sim;
          SH_EXP2_SUM   [ thread_id] += tmp_exp * tmp_exp;
          SH_SIM2_SUM   [ thread_id] += tmp_sim * tmp_sim;
          SH_COUNT_SUM  [ thread_id] += 1;
        }

        // go to next block
        index_z += grid_size_z;
      }
      // go to next block
      index_y += grid_size_y;
    }
    // go to next block
    index_x += grid_size_x;
  }

  // wait for block
  barrier( CLK_LOCAL_MEM_FENCE);

  // do reduction in shared memory for each group
  // reduce within block
  for( size_t offset = block_size / 2; offset > 0; offset = offset / 2)
  {
    if( thread_id < offset)
    {
      SH_EXP_SUM    [ thread_id] += SH_EXP_SUM    [ thread_id + offset];
      SH_SIM_SUM    [ thread_id] += SH_SIM_SUM    [ thread_id + offset];
      SH_EXP_SIM_SUM[ thread_id] += SH_EXP_SIM_SUM[ thread_id + offset];
      SH_EXP2_SUM   [ thread_id] += SH_EXP2_SUM   [ thread_id + offset];
      SH_SIM2_SUM   [ thread_id] += SH_SIM2_SUM   [ thread_id + offset];
      SH_COUNT_SUM  [ thread_id] += SH_COUNT_SUM  [ thread_id + offset];
      barrier( CLK_LOCAL_MEM_FENCE);
    }
  }

  // copy result to output
  if( thread_id == 0)
  {
    EXP_SUM_OUTPUT    [ get_group_id( 2) * get_num_groups( 1) * get_num_groups( 0) + get_group_id( 1) * get_num_groups( 0) + get_group_id( 0)] = SH_EXP_SUM    [ 0];
    SIM_SUM_OUTPUT    [ get_group_id( 2) * get_num_groups( 1) * get_num_groups( 0) + get_group_id( 1) * get_num_groups( 0) + get_group_id( 0)] = SH_SIM_SUM    [ 0];
    EXP_SIM_SUM_OUTPUT[ get_group_id( 2) * get_num_groups( 1) * get_num_groups( 0) + get_group_id( 1) * get_num_groups( 0) + get_group_id( 0)] = SH_EXP_SIM_SUM[ 0];
    EXP2_SUM_OUTPUT   [ get_group_id( 2) * get_num_groups( 1) * get_num_groups( 0) + get_group_id( 1) * get_num_groups( 0) + get_group_id( 0)] = SH_EXP2_SUM   [ 0];
    SIM2_SUM_OUTPUT   [ get_group_id( 2) * get_num_groups( 1) * get_num_groups( 0) + get_group_id( 1) * get_num_groups( 0) + get_group_id( 0)] = SH_SIM2_SUM   [ 0];
    COUNT_OUTPUT      [ get_group_id( 2) * get_num_groups( 1) * get_num_groups( 0) + get_group_id( 1) * get_num_groups( 0) + get_group_id( 0)] = SH_COUNT_SUM  [ 0];
  }
}
