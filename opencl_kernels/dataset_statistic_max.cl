// this is only included, so that eclipse c++ editor can be used and parses information correctly
// this will not be interpreted, when the program is compiled with an OpenCL compiler
#ifndef __OPENCL_VERSION__
 #include "kernel.h"
#endif

//! @brief determine the maximal value of each vector elements in a data set
//! @param DATASET the dataset with NUM_DATA_ITEMS vectors and VECTOR_LENGTH elements each
//! @param RESULT_MATRIX the matrix that will contain the unreduced results for each column in a row
//! @param NUM_DATA_ITEMS number of rows in dataset
//! @param VECTOR_LENGTH length of each vector
//! @param SHARED shared memory for the reduction
__kernel
void DataSetStatisticMax
(
  __const __global PRECISION* DATASET,
  __const          uint       NUM_DATA_ITEMS,
  __const          uint       VECTOR_LENGTH,
          __global PRECISION* RESULT_MATRIX,
          __local  PRECISION* SHARED
)
{
  __const size_t vector_index = get_global_id( 1);
  size_t global_index = get_global_id( 0);

  PRECISION accumulator = -INFINITY;

  // Loop sequentially over chunks of input vector
  while( global_index < NUM_DATA_ITEMS)
  {
    accumulator = max( accumulator, DATASET[ global_index * VECTOR_LENGTH + vector_index]);
    global_index += get_global_size( 0);
  }

  // Perform parallel reduction
  __const size_t local_index = get_local_id( 0);
  SHARED[ local_index] = accumulator;
  barrier( CLK_LOCAL_MEM_FENCE);

  for( size_t offset = get_local_size( 0) / 2; offset > 0; offset = offset / 2)
  {
    if( local_index < offset)
    {
      SHARED[ local_index] = max( SHARED[ local_index + offset], SHARED[ local_index]);
    }
    barrier( CLK_LOCAL_MEM_FENCE);
  }
  if( local_index == 0)
  {
    RESULT_MATRIX[ vector_index * get_num_groups( 0) + get_group_id( 0)] = SHARED[ 0];
  }
}
