// this is only included, so that eclipse c++ editor can be used and parses information correctly
// this will not be interpreted, when the program is compiled with an OpenCL compiler
#ifndef __OPENCL_VERSION__
 #include "kernel.h"
#endif

//! @brief calculates partial rmsd as the resulting array is num_groups size
//! @param INPUT_A      calculated
//! @param INPUT_B      actual
//! @param OUTPUT       output of temp reduction array of size number groups
//! @param NUM_ELEMENTS number of elements in the array
//! @param SHARED       __local mem allocation of size block_size * sizeof( precision)
__kernel void RmsdKernel
(
  const __global PRECISION* INPUT_A,           // calculated
  const __global PRECISION* INPUT_B,           // actual
        __global PRECISION* OUTPUT,            // output of temp reduction array of size number groups
  const          int        NUM_ELEMENTS,      // number of elements in the array
        __local  PRECISION* SHARED             // __local mem allocation of size block_size * sizeof( precision)
)
{
  // acquire position in grid and sizes
  const unsigned int block_size = get_local_size( 0);
  const unsigned int thread_id  = get_local_id( 0);
  const unsigned int grid_size  = block_size * get_num_groups( 0);

  // start index for that thread
  unsigned int index            = get_group_id( 0) * block_size + get_local_id( 0);

  // temporary variable that holds value from global memory
  PRECISION tmp_a = 0;
  PRECISION tmp_b = 0;
  PRECISION diff  = 0;

  // Initialize elements in shared memory to 0
  SHARED[ thread_id] = 0;

  // perform first level of reduction,
  // reading from global memory, writing to shared memory
  while( index < NUM_ELEMENTS)
  {
    // simulated intensity
    tmp_a = INPUT_A[ index];
    tmp_b = INPUT_B[ index];
    diff = tmp_a - tmp_b;
    SHARED[ thread_id] += diff * diff;

    // go to next block
    index += grid_size;
  }

  // wait for block
  barrier( CLK_LOCAL_MEM_FENCE);

  // do reduction in shared memory
  // reduce within block
  for( int offset = block_size / 2; offset > 0; offset = offset / 2)
  {
    if( thread_id < offset)
    {
      SHARED[ thread_id] += SHARED[ thread_id + offset];
    }
    barrier( CLK_LOCAL_MEM_FENCE);
  }

  // copy result to output
  if( thread_id == 0)
  {
    OUTPUT[ get_group_id( 0)] = SHARED[ 0];
  }
}
