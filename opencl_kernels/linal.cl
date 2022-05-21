// this is only included, so that eclipse c++ editor can be used and parses information correctly
// this will not be interpreted, when the program is compiled with an OpenCL compiler
#ifndef __OPENCL_VERSION__
 #include "kernel.h"
#endif

#pragma OPENCL FP_CONTRACT OFF

//! @brief multiply two matrices
//! @param MATRIX_A first matrix
//! @param MATRIX_B second matrix
//! @param NUM_COLS_A number of columns in MATRIX_A
//! @param NUM_COLS_B number of columns in MATRIX_B
//! @param OUTPUT_MATRIX result matrix
//! @param SHARED_A shared memory for tile of MATRIX_A
//! @param SHARED_b shared memory for tile of MATRIX_B
__kernel void MatrixMultiplication
(
  const __global PRECISION* MATRIX_A,
  const __global PRECISION* MATRIX_B,
  const          uint       NUM_COLS_A,
  const          uint       NUM_COLS_B,
        __global PRECISION* OUTPUT_MATRIX,
        __local  PRECISION* SHARED_A,
        __local  PRECISION* SHARED_B
)
{
  // Block index
  const size_t bx = get_group_id( 0);
  const size_t by = get_group_id( 1);

  // Thread index
  const size_t tx = get_local_id( 0);
  const size_t ty = get_local_id( 1);

  const size_t block_size = get_local_size( 0);

  // Index of the first sub-matrix of MATRIX_A processed by the block
  const size_t a_beg = NUM_COLS_A * block_size * by;

  // Index of the end sub-matrix of MATRIX_A processed by the block
  const size_t a_end = a_beg + NUM_COLS_A - 1;

  // Step size used to iterate through the sub-matrices of MATRIX_A
  const size_t a_step = block_size;

  // Index of the first sub-matrix of MATRIX_B processed by the block
  const size_t b_beg = block_size * bx;

  // Step size used to iterate through the sub-matrices of MATRIX_B
  const size_t b_step = block_size * NUM_COLS_B;

  // Csub is used to store the element of the block sub-matrix
  // that is computed by the thread
  PRECISION c_sub = 0.0;

  // Loop over all the sub-matrices of MATRIX_A and MATRIX_B
  // required to compute the block sub-matrix
  for( size_t a = a_beg, b = b_beg; a <  a_end; a += a_step, b += b_step)
  {
#define AS( i, j) SHARED_A[ i * block_size + j]
#define BS( i, j) SHARED_B[ i * block_size + j]
    // Load the matrices from device memory
    // to shared memory; each thread loads
    // one element of each matrix
    AS( ty, tx) = MATRIX_A[ a + NUM_COLS_A * ty + tx];
    BS( ty, tx) = MATRIX_B[ b + NUM_COLS_B * ty + tx];
    // Synchronize to make sure the matrices are loaded
    barrier( CLK_LOCAL_MEM_FENCE);

    // Multiply the two matrices together;
    // each thread computes one element
    // of the block sub-matrix
//    #pragma unroll
    for( size_t k = 0; k < block_size; ++k)
    {
      c_sub += AS( ty, k) * BS( k, tx);
    }
#undef AS
#undef BS
    // Synchronize to make sure that the preceding
    // computation is done before loading two new
    // sub-matrices of MATRIX_A and MATRIX_B in the next iteration
    barrier( CLK_LOCAL_MEM_FENCE);
  }

  // Write the block sub-matrix to device memory;
  // each thread writes one element
  OUTPUT_MATRIX[ get_global_id( 1) * get_global_size( 0) + get_global_id( 0)] = c_sub;
}

//! @brief multiply two matrices
//! @param MATRIX_A first matrix
//! @param MATRIX_B second matrix
//! @param NUM_COLS_A number of columns in MATRIX_A
//! @param NUM_COLS_B number of columns in MATRIX_B
//! @param OUTPUT_MATRIX result matrix
//! @param SHARED_A shared memory for tile of MATRIX_A
//! @param SHARED_b shared memory for tile of MATRIX_B
__kernel void MatrixMultiplicationUnrolled
(
  const __global PRECISION* MATRIX_A,
  const __global PRECISION* MATRIX_B,
  const          uint       NUM_COLS_A,
  const          uint       NUM_COLS_B,
        __global PRECISION* OUTPUT_MATRIX,
        __local  PRECISION* SHARED_A,
        __local  PRECISION* SHARED_B
)
{
  // Block index
  const size_t bx = get_group_id( 0);
  const size_t by = get_group_id( 1);

  // Thread index
  const size_t tx = get_local_id( 0);
  const size_t ty = get_local_id( 1);

  const size_t block_size = get_local_size( 0);

  // Index of the first sub-matrix of MATRIX_A processed by the block
  const size_t a_beg = NUM_COLS_A * block_size * by;

  // Index of the end sub-matrix of MATRIX_A processed by the block
  const size_t a_end = a_beg + NUM_COLS_A - 1;

  // Step size used to iterate through the sub-matrices of MATRIX_A
  const size_t a_step = block_size;

  // Index of the first sub-matrix of MATRIX_B processed by the block
  const size_t b_beg = block_size * bx;

  // Step size used to iterate through the sub-matrices of MATRIX_B
  const size_t b_step = block_size * NUM_COLS_B;

  // Csub is used to store the element of the block sub-matrix
  // that is computed by the thread
  float4 c_sub = ( float4)( 0.0, 0.0, 0.0, 0.0);

  // Loop over all the sub-matrices of MATRIX_A and MATRIX_B
  // required to compute the block sub-matrix
  for( size_t a = a_beg, b = b_beg; a <  a_end; a += a_step, b += b_step)
  {
#define AS( i, j) SHARED_A[ i * block_size + j]
#define BS( i, j) SHARED_B[ i * block_size + j]
    // Load the matrices from device memory
    // to shared memory; each thread loads
    // one element of each matrix
    AS( ty, tx) = MATRIX_A[ a + NUM_COLS_A * ty + tx];
    BS( ty, tx) = MATRIX_B[ b + NUM_COLS_B * ty + tx];
    // Synchronize to make sure the matrices are loaded
    barrier( CLK_LOCAL_MEM_FENCE);

    // Multiply the two matrices together;
    // each thread computes one element
    // of the block sub-matrix
    float4 temp_a0 = ( float4)( AS( ty, 0), AS( ty, 1), AS( ty, 2), AS( ty, 3));
    float4 temp_b0 = ( float4)( BS( 0, tx), BS( 1, tx), BS( 2, tx), BS( 3, tx));
    c_sub = mad( temp_a0, temp_b0, c_sub);
    float4 temp_a1 = ( float4)( AS( ty, 4), AS( ty, 5), AS( ty, 6), AS( ty, 7));
    float4 temp_b1 = ( float4)( BS( 4, tx), BS( 5, tx), BS( 6, tx), BS( 7, tx));
    c_sub = mad( temp_a1, temp_b1, c_sub);
    float4 temp_a2 = ( float4)( AS( ty, 8), AS( ty, 9), AS( ty, 10), AS( ty, 11));
    float4 temp_b2 = ( float4)( BS( 8, tx), BS( 9, tx), BS( 10, tx), BS( 11, tx));
    c_sub = mad( temp_a2, temp_b2, c_sub);
    float4 temp_a3 = ( float4)( AS( ty, 12), AS( ty, 13), AS( ty, 14), AS( ty, 15));
    float4 temp_b3 = ( float4)( BS( 12, tx), BS( 13, tx), BS( 14, tx), BS( 15, tx));
    c_sub = mad( temp_a3, temp_b3, c_sub);
#undef AS
#undef BS

    // Synchronize to make sure that the preceding
    // computation is done before loading two new
    // sub-matrices of MATRIX_A and MATRIX_B in the next iteration
    barrier( CLK_LOCAL_MEM_FENCE);
  }
  // Write the block sub-matrix to device memory;
  // each thread writes one element
  OUTPUT_MATRIX[ get_global_id( 1) * get_global_size( 0) + get_global_id( 0)] = c_sub.x + c_sub.y + c_sub.z + c_sub.w;
}
#pragma OPENCL FP_CONTRACT ON

//! @brief multiply vector with MATRIX resulting in a vector
//! @param MATRIX row major matrix to multiply
//! @param VECTOR vector to multiply
//! @param WIDTH width of matrix
//! @param HEIGHT height of matrix
//! @param RESULT result vector
//! @param SHARED shared memory
__kernel void MatrixVectorMultiplication
(
  const __global PRECISION* MATRIX,
  const __global PRECISION* VECTOR,
  const          uint       WIDTH,
  const          uint       HEIGHT,
        __global PRECISION* RESULT,
        __local  PRECISION* SHARED
)
{
  // Each work-group computes multiple elements of RESULT
  for( size_t y = get_group_id( 0); y < HEIGHT; y += get_num_groups( 0))
  {
    const __global PRECISION* row = MATRIX + y * WIDTH;

    // Each work-item accumulates as many products as necessary
    // into local variable sum
    PRECISION sum = 0;
    for( size_t x = get_local_id( 0); x < WIDTH; x += get_local_size( 0))
    {
      sum += row[ x] * VECTOR[ x];
    }

    // Each partial dot product is stored in shared memory
    SHARED[ get_local_id( 0)] = sum;

    // Perform parallel reduction to add each work-item's
    // partial dot product together
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
        // Add two elements from the SHARED array
        // and store the result in SHARED[index]
        SHARED[ get_local_id( 0)] += SHARED[ get_local_id( 0) + stride];
      }
    }
    // Write the result of the reduction to global memory
    if( get_local_id( 0) == 0)
    {
      RESULT[ y] = SHARED[ 0];
    }

    // Synchronize to make sure each warp is done reading
    // SHARED before it is overwritten in the next step
    barrier( CLK_LOCAL_MEM_FENCE);
  }
}

//! @brief multiply vector with MATRIX resulting in a vector
//! @param MATRIX row major matrix to multiply
//! @param VECTOR vector to multiply
//! @param WIDTH width of matrix
//! @param HEIGHT height of matrix
//! @param RESULT result vector
//! @param SHARED shared memory
__kernel void MatrixVectorMultiplicationPlusVector
(
  const __global PRECISION* MATRIX,
  const __global PRECISION* VECTOR,
  const __global PRECISION* ADD_VECTOR,
  const          uint       WIDTH,
  const          uint       HEIGHT,
        __global PRECISION* RESULT,
        __local  PRECISION* SHARED
)
{
  // Each work-group computes multiple elements of RESULT
  for( size_t y = get_group_id( 0); y < HEIGHT; y += get_num_groups( 0))
  {
    const __global PRECISION* row = MATRIX + y * WIDTH;

    // Each work-item accumulates as many products as necessary
    // into local variable sum
    PRECISION sum = 0;
    for( size_t x = get_local_id( 0); x < WIDTH; x += get_local_size( 0))
    {
      sum += row[ x] * VECTOR[ x];
    }

    // Each partial dot product is stored in shared memory
    SHARED[ get_local_id( 0)] = sum;

    // Perform parallel reduction to add each work-item's
    // partial dot product together
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
        // Add two elements from the SHARED array
        // and store the result in SHARED[index]
        SHARED[ get_local_id( 0)] += SHARED[ get_local_id( 0) + stride];
      }
    }
    // Write the result of the reduction to global memory
    if( get_local_id( 0) == 0)
    {
      RESULT[ y] = SHARED[ 0] + ADD_VECTOR[ y];
    }

    // Synchronize to make sure each warp is done reading
    // SHARED before it is overwritten in the next step
    barrier( CLK_LOCAL_MEM_FENCE);
  }
}

//! @brief dot product of two vectors
//! @param VECTOR_A first vector
//! @param VECTOR_B second vector
//! @param LENGTH number of elements in vectors to consider
//! @param RESULT the vector containing the element-wise products
//! @param SHARED shared memory for reduction
__kernel void VectorDotProduct
(
  const __global PRECISION* VECTOR_A,
  const __global PRECISION* VECTOR_B,
  const          uint       LENGTH,
        __global PRECISION* RESULT,
        __local  PRECISION* SHARED
)
{
  // perform first level of reduction,
  // reading from global memory, writing to shared memory
  const uint thread_id = get_local_id( 0);
  const uint block_size = get_local_size( 0);
  const uint grid_size = block_size * 2 * get_num_groups( 0);
  uint index = get_group_id( 0) * block_size * 2 + thread_id;
  SHARED[ thread_id] = 0;

  // we reduce multiple elements per thread.  The number is determined by the
  // number of active thread blocks (via gridDim).  More blocks will result
  // in a larger grid_size and therefore fewer elements per thread
  while( index < LENGTH)
  {
    SHARED[ thread_id] += VECTOR_A[ index] * VECTOR_B[ index];
    // ensure we don't read out of bounds
    if( index + block_size < LENGTH)
    {
      SHARED[ thread_id] += VECTOR_A[ index + block_size] * VECTOR_B[ index + block_size];
    }
    index += grid_size;
  }

  barrier( CLK_LOCAL_MEM_FENCE);

  // do reduction in shared mem
  for( uint offset = block_size / 2; offset > 0; offset = offset / 2)
  {
    if( thread_id < offset)
    {
      SHARED[ thread_id] += SHARED[ thread_id + offset];
    }
    barrier( CLK_LOCAL_MEM_FENCE);
  }

  // write result for this block to global mem
  if( thread_id == 0)
  {
    RESULT[ get_group_id( 0)] = SHARED[ 0];
  }
  barrier( CLK_LOCAL_MEM_FENCE);
}

//! @brief inner product of a vector
//! @param VECTOR the vector
//! @param LENGTH number of elements in vector
//! @param RESULT vector of inner products
//! @param SHARED shared memory for reduction
__kernel void VectorInnerProduct
(
  const __global PRECISION* VECTOR,
  const          uint       LENGTH,
        __global PRECISION* RESULT,
        __local  PRECISION* SHARED
)
{
  uint global_index = get_global_id( 0);
  uint greater_than = 0;
  PRECISION accumulator = 0;

  // Loop sequentially over chunks of input vector
  while( global_index < LENGTH)
  {
    PRECISION element = VECTOR[ global_index];
    accumulator += element * element;
    global_index += get_global_size( 0);
  }

  // Perform parallel reduction
  uint local_index = get_local_id( 0);
  SHARED[ local_index] = accumulator;
  barrier( CLK_LOCAL_MEM_FENCE);

  for( int offset = get_local_size( 0) / 2; offset > 0; offset = offset / 2)
  {
    if( local_index < offset)
    {
      SHARED[ local_index] += SHARED[ local_index + offset];
    }
    barrier( CLK_LOCAL_MEM_FENCE);
  }
  if( local_index == 0)
  {
    RESULT[ get_group_id( 0)] = SHARED[ 0];
  }
}

//! @brief sum of elements in a vector
//! @param INPUT input vector
//! @param NUM_ELEMENTS number of elements in VECTOR
//! @param OUTPUT vector with partial sums of input vector, still needs to be reduced
//! @param SHARED shared memory for reduction
__kernel void ReductionSum
(
  const __global PRECISION* INPUT,
  const          uint       NUM_ELEMENTS,
        __global PRECISION* OUTPUT,
        __local  PRECISION* SHARED
)
{
  // perform first level of reduction,
  // reading from global memory, writing to shared memory
  const uint thread_id = get_local_id( 0);
  const uint block_size = get_local_size( 0);
  const uint grid_size = block_size * 2 * get_num_groups( 0);
  uint index = get_group_id( 0) * get_local_size( 0) * 2 + get_local_id( 0);
  SHARED[ thread_id] = 0;

  // we reduce multiple elements per thread.  The number is determined by the
  // number of active thread blocks (via gridDim).  More blocks will result
  // in a larger grid_size and therefore fewer elements per thread
  while( index < NUM_ELEMENTS)
  {
    SHARED[ thread_id] += INPUT[ index];
    // ensure we don't read out of bounds
    if( index + block_size < NUM_ELEMENTS)
    {
      SHARED[ thread_id] += INPUT[ index + block_size];
    }
    index += grid_size;
  }

  barrier( CLK_LOCAL_MEM_FENCE);

  // do reduction in shared mem
  for( uint offset = block_size / 2; offset > 0; offset = offset / 2)
  {
    if( thread_id < offset)
    {
      SHARED[ thread_id] += SHARED[ thread_id + offset];
    }
    barrier( CLK_LOCAL_MEM_FENCE);
  }

  // write result for this block to global mem
  if( thread_id == 0)
  {
    OUTPUT[ get_group_id( 0)] = SHARED[ 0];
  }
}

//! @brief transposes input matrix
//! @param INPUT  input
//! @param ROWS   number rows
//! @param COLS   number cols
//! @param OUTPUT transposed output
//! @param SHARED  local mem allocation
__kernel void MatrixTranspose
(
  const __global PRECISION *INPUT,   // input
  const          uint       ROWS,    // number rows
  const          uint       COLS,    // number cols
        __global PRECISION *OUTPUT,  // transposed output
        __local  PRECISION *SHARED   // local mem allocation
)
{
  const int block_size = get_local_size( 0);

  // read the matrix tile into shared memory
  uint index_x = get_global_id( 0);
  uint index_y = get_global_id( 1);

  if(( index_x < COLS) && ( index_y < ROWS))
  {
    const uint index_in = index_y * COLS + index_x;
    SHARED[ get_local_id( 1) * block_size + get_local_id( 0)] = INPUT[ index_in];
  }
  barrier( CLK_LOCAL_MEM_FENCE);

  // write the transposed matrix tile to global memory
  index_x = get_group_id( 1) * block_size + get_local_id( 0);
  index_y = get_group_id( 0) * block_size + get_local_id( 1);
  if(( index_x < ROWS) && ( index_y < COLS))
  {
    const uint index_out = index_y * ROWS + index_x;
    OUTPUT[ index_out] = SHARED[ get_local_id( 0) * block_size + get_local_id( 1)];
  }
}

__kernel void ArgMax
(
  const __global PRECISION* INPUT,
        __global int*       INDEX,
        __global PRECISION* RESULT,
  const          uint       NUM_ELEMENTS,
        __local  PRECISION* SHARED,
        __local  int*       IND_SHARED
)
{
  uint global_index = get_global_id( 0);
  uint greater_than = 0;
  PRECISION accumulator = -INFINITY;

  // Loop sequentially over chunks of input vector
  while( global_index < NUM_ELEMENTS)
  {
    PRECISION element = INPUT[ global_index];
    greater_than = isgreater( accumulator, element);
    accumulator = greater_than ? accumulator : element;
    global_index += get_global_size( 0);
  }

  // Perform parallel reduction
  uint local_index = get_local_id( 0);
  SHARED[ local_index] = accumulator;
  IND_SHARED[ local_index] = get_global_id( 0);
  barrier( CLK_LOCAL_MEM_FENCE);

  for( int offset = get_local_size( 0) / 2; offset > 0; offset = offset / 2)
  {
    if( local_index < offset)
    {
      greater_than = isgreater( SHARED[ local_index + offset], SHARED[ local_index]);
      SHARED[ local_index] = greater_than ? SHARED[ local_index + offset] : SHARED[ local_index];
      IND_SHARED[ local_index] = greater_than ? IND_SHARED[ local_index + offset] : IND_SHARED[ local_index];
    }
    barrier( CLK_LOCAL_MEM_FENCE);
  }
  if( local_index == 0)
  {
    RESULT[ get_group_id( 0)] = SHARED[ 0];
    INDEX[ get_group_id( 0)] = IND_SHARED[ 0];
  }
}

__kernel void ArgMin
(
  const __global PRECISION* INPUT,
        __global int*       INDEX,
        __global PRECISION* RESULT,
  const          uint       NUM_ELEMENTS,
        __local  PRECISION* SHARED,
        __local  int*       IND_SHARED
)
{
  int global_index = get_global_id( 0);
  int less_than = 0;
  int index = 0;
  PRECISION accumulator = INFINITY;

  // Loop sequentially over chunks of input vector
  while( global_index < NUM_ELEMENTS)
  {
    PRECISION element = INPUT[ global_index];
    less_than = isless( accumulator, element);
    accumulator = less_than ? accumulator : element;
    index = less_than ? 0 : global_index;
    global_index += get_global_size( 0);
  }

  // Perform parallel reduction
  int local_index = get_local_id( 0);
  SHARED[ local_index] = accumulator;
  IND_SHARED[ local_index] = get_global_id( 0);
  barrier( CLK_LOCAL_MEM_FENCE);

  for( int offset = get_local_size( 0) / 2; offset > 0; offset = offset / 2)
  {
    if( local_index < offset)
    {
      less_than = isless( SHARED[ local_index + offset], SHARED[ local_index]);
      SHARED[ local_index] = less_than ? SHARED[ local_index + offset] : SHARED[ local_index];
      IND_SHARED[ local_index] = less_than ? IND_SHARED[ local_index + offset] : IND_SHARED[ local_index];
    }
    barrier( CLK_LOCAL_MEM_FENCE);
  }
  if( local_index == 0)
  {
    RESULT[ get_group_id( 0)] = SHARED[ 0];
    INDEX[ get_group_id( 0)] = IND_SHARED[ 0];
  }
}

__kernel void ReductionMax
(
  const __global PRECISION* INPUT,
        __global PRECISION* RESULT,
  const          uint       NUM_ELEMENTS,
        __local  PRECISION* SHARED
)
{
  uint global_index = get_global_id( 0);
  PRECISION accumulator = 0;

  // Loop sequentially over chunks of input vector
  while( global_index < NUM_ELEMENTS)
  {
    PRECISION element = INPUT[ global_index];
    accumulator = ( accumulator > element) ? accumulator : element;
    global_index += get_global_size( 0);
  }

  // Perform parallel reduction
  uint local_index = get_local_id( 0);
  SHARED[ local_index] = accumulator;
  barrier( CLK_LOCAL_MEM_FENCE);
  for( int offset = get_local_size( 0) / 2; offset > 0; offset = offset / 2)
  {
    if( local_index < offset)
    {
      PRECISION other = SHARED[ local_index + offset];
      PRECISION mine = SHARED[ local_index];
      SHARED[ local_index] = ( mine > other) ? mine : other;
    }
    barrier( CLK_LOCAL_MEM_FENCE);
  }
  if( local_index == 0)
  {
    RESULT[ get_group_id( 0)] = SHARED[ 0];
  }
}

__kernel void ReductionMin
(
  const __global PRECISION* INPUT,
        __global PRECISION* RESULT,
  const          uint       NUM_ELEMENTS,
        __local  PRECISION* SHARED
)
{
  uint global_index = get_global_id( 0);
  PRECISION accumulator = INFINITY;
  // Loop sequentially over chunks of input vector
  while( global_index < NUM_ELEMENTS)
  {
    PRECISION element = INPUT[ global_index];
    accumulator = ( accumulator < element) ? accumulator : element;
    global_index += get_global_size( 0);
  }

  // Perform parallel reduction
  int local_index = get_local_id( 0);
  SHARED[ local_index] = accumulator;
  barrier( CLK_LOCAL_MEM_FENCE);
  for( int offset = get_local_size( 0) / 2; offset > 0; offset = offset / 2)
  {
    if( local_index < offset)
    {
      PRECISION other = SHARED[ local_index + offset];
      PRECISION mine = SHARED[ local_index];
      SHARED[ local_index] = ( mine < other) ? mine : other;
    }
    barrier( CLK_LOCAL_MEM_FENCE);
  }
  if( local_index == 0)
  {
    RESULT[ get_group_id( 0)] = SHARED[ 0];
  }
}


//! @brief  elementwise addition kernel
//! @param  MATRIX_ADDED_TO the weights which will be added to
//! @param  MATRIX_TO_ADD the values to be added to the weights
//! @param  ROWS   rows in input matrix
//! @param  COLS   cols in input matrix
//! @param  XPAD   padding added to cols
//! @param  YPAD   padding added to rows
__kernel void MatrixAddition
(
        __global PRECISION *MATRIX_ADDED_TO,
  const __global PRECISION *MATRIX_TO_ADD,
  const          uint       ROWS,
  const          uint       COLS,
  const          uint       XPAD,
  const          uint       YPAD
)
{
  const uint index_x = get_global_id( 0);
  const uint index_y = get_global_id( 1);
  const uint size_x = COLS - XPAD;
  const uint size_y = ROWS - YPAD;
  const uint index = index_y * COLS + index_x;
  if( index_x < size_x && index_y < size_y)
  {
    MATRIX_ADDED_TO[ index] += MATRIX_TO_ADD[ index];
  }
}

//! @brief  adds vector to each row in a matrix
//! @param  VECTOR  the vector to be added to the input rows
//! @param  MATRIX  the input matrix
//! @param  ROWS   rows in input matrix
//! @param  COLS   cols in input matrix
//! @param  XPAD   padding added to cols
//! @param  YPAD   padding added to rows
__kernel void VectorMatrixAddition
(
  const __global PRECISION *VECTOR,
        __global PRECISION *MATRIX,
  const          uint       ROWS,
  const          uint       COLS,
  const          uint       XPAD,
  const          uint       YPAD
)
{
  const uint index_x = get_global_id( 0);
  const uint index_y = get_global_id( 1);
  const uint boundary_x = COLS - XPAD;
  const uint boundary_y = ROWS - YPAD;
  const uint index = index_y * COLS + index_x;
  if( index_x >= boundary_x || index_y >= boundary_y)
  {
    return;
  }

  MATRIX[ index] += VECTOR[ index_x];
}

__kernel void InvertDiagonalMatrix
(
        __global PRECISION *INPUT,
  const          uint       ROWS,
  const          uint       ROW_PAD
)
{
  const uint i = get_global_id( 0);
  if( i >= ROWS) return;

  if( INPUT[ i * ( ROWS + ROW_PAD) + i] == ( PRECISION)( 0))
  {
    return;
  }

  INPUT[ i * ( ROWS + ROW_PAD) + i] = ( PRECISION)( 1) / INPUT[ i * ( ROWS + ROW_PAD) + i];
}

//! @brief determinant of a 3x3 matrix
//! @param MATRIX the matrix of local_size( 0) x localsize( 1)
//! @param DETERMINANT the result
__kernel void Matrix3x3Determinant
(
  const __global PRECISION *MATRIX,
        __global PRECISION *DETERMINANT,
        __local  PRECISION *SHARED_MATRIX,
        __local  PRECISION *SHARED_RESULT
)
{
  const uint thread_x = get_local_id( 0);
  const uint thread_y = get_local_id( 1);
  const uint block_size = get_local_size( 0);

  SHARED_MATRIX[ thread_y * block_size + thread_x] = MATRIX[ thread_y * block_size + thread_x];
  if( thread_y != 0)
  {
    return;
  }

  SHARED_RESULT[ thread_x] = 0.0;
  // synchronize
  barrier( CLK_LOCAL_MEM_FENCE);

#define matrix( i, j) SHARED_MATRIX[ i * block_size + j]

  const uint length = 4;
  if( thread_x < 3)
  {
    SHARED_RESULT[ thread_x         ] =  matrix( 0,   thread_x         ) * matrix( 1, ( thread_x + 1) % 3) * matrix( 2, ( thread_x + 2) % 3);
    SHARED_RESULT[ thread_x + length] = -matrix( 0, ( thread_x + 2) % 3) * matrix( 1, ( thread_x + 1) % 3) * matrix( 2, ( thread_x + 0) % 3);
  }
  barrier( CLK_LOCAL_MEM_FENCE);

  // reduction
  for( uint offset = length; offset > 0; offset = offset / 2)
  {
    if( thread_x < offset)
    {
      SHARED_RESULT[ thread_x] += SHARED_RESULT[ thread_x + offset];
    }
    barrier( CLK_LOCAL_MEM_FENCE);
  }
  if( thread_x == 0)
  {
    DETERMINANT[ 0] = SHARED_RESULT[ 0];
  }
#undef matrix
}

//! @brief eigenvalues of a symmetric matrix with 3x3 elements in use
//! @param MATRIX the matrix of size localsize( 0) x localsize( 1)
//! @param SQRT take the squareroot of each eigenvalue, which can be used (true) if matrix was symmetrized before
//! @param EIGENVALUES the result eigenvalues
void Matrix3x3SymmetricEigenvaluesAux
(
  const PRECISION *MATRIX,
  const      uint  SQRT,
        PRECISION *EIGENVALUES
)
{
#define matrix( i, j) MATRIX[ i * 4 + j]
  // Determine coefficients of characteristic poynomial. We write
  //       | a   d   f  |
  //  A =  | d*  b   e  |
  //       | f*  e*  c  |
  const PRECISION de = matrix( 0, 1) * matrix( 1, 2); // d * e
  const PRECISION dd = pown( matrix( 0, 1), 2);       // d^2
  const PRECISION ee = pown( matrix( 1, 2), 2);       // e^2
  const PRECISION ff = pown( matrix( 0, 2), 2);       // f^2
  const PRECISION  m = matrix( 0, 0) + matrix( 1, 1) + matrix( 2, 2);
  const PRECISION c1 = matrix( 0, 0) * matrix( 1, 1) + matrix( 0, 0) * matrix( 2, 2) + matrix( 1, 1) * matrix( 2, 2) - dd - ee - ff; // a*b + a*c + b*c - d^2 - e^2 - f^2
  const PRECISION c0 = matrix( 2, 2) * dd            + matrix( 0, 0) * ee            + matrix( 1, 1) * ff            - matrix( 0, 0) * matrix( 1, 1) * matrix( 2, 2) - 2 * matrix( 0, 2) * de; // c*d^2 + a*e^2 + b*f^2 - a*b*c - 2*f*d*e)

  const PRECISION p = pown( m, 2) - 3 * c1;
  const PRECISION q = m * ( p - ( 3.0 / 2.0) * c1) - ( 27.0 / 2.0) * c0;
  const PRECISION sqrt_p = sqrt( fabs( p));

  PRECISION phi = 27.0 * ( 0.25 * pown( c1, 2) * ( p - c1) + c0 * ( q + 27.0 / 4.0 * c0));
  phi = ( 1.0 / 3.0) * atan2( sqrt( fabs( phi)), q);

  const PRECISION c = sqrt_p * cos( phi);
  const PRECISION s = ( 1.0 / sqrt( 3.0)) * sqrt_p * sin( phi);

  const PRECISION tmp = ( 1.0 / 3.0) * ( m - c);

  if( SQRT == 0)
  {
    EIGENVALUES[ 0] = tmp + c;
    EIGENVALUES[ 1] = tmp - s;
    EIGENVALUES[ 2] = tmp + s;
  }
  else
  {
    EIGENVALUES[ 0] = sqrt( tmp + c);
    EIGENVALUES[ 1] = sqrt( tmp - s);
    EIGENVALUES[ 2] = sqrt( tmp + s);
  }
#undef matrix
}

//! @brief eigenvalues of a symmetric matrix with 3x3 elements in use
//! @param MATRIX the matrix of size localsize( 0) x localsize( 1)
//! @param SQRT take the squareroot of each eigenvalue, which can be used (true) if matrix was symmetrized before
//! @param EIGENVALUES the result eigenvalues
__kernel void Matrix3x3SymmetricEigenvalues
(
  const __global PRECISION *MATRIX,
  const               uint  SQRT,
        __global PRECISION *EIGEN_VALUES
)
{
  PRECISION matrix_local[ 16];
  PRECISION eigen_values_local[ 4];

  // copy input matrix
  for( uint row = 0; row < 4; ++row)
  {
    eigen_values_local[ row] = 0.0;
    for( uint col = 0; col < 4; ++col)
    {
      matrix_local[ row * 4 + col] = MATRIX[ row * 16 + col];
    }
  }

  // calculate with auxilary functions
  Matrix3x3SymmetricEigenvaluesAux( matrix_local, SQRT, eigen_values_local);

  // copy results
  for( uint row = 0; row < 4; ++row)
  {
    EIGEN_VALUES[ row] = eigen_values_local[ row];
  }
}

//! @brief transform matrix to tridiagonal form
// Reduces a symmetric 3x3 matrix to tridiagonal form by applying
// (unitary) Householder transformations:
//            [ d[0]  e[0]       ]
//    A = Q . [ e[0]  d[1]  e[1] ] . Q^T
//            [       e[1]  d[2] ]
//! @param EIGEN_VECTORS storage for Q
//! @param EIGENVALUES storage for diagonal elements d
//! @param OFF_DIAGONAL storage for off diagonal elements
void Matrix3x3SymmetricTridiagonalizeHouseholderAux
(
  const PRECISION *MATRIX,
        PRECISION *EIGEN_VECTORS,
        PRECISION *EIGENVALUES,
        PRECISION *OFF_DIAGONAL
)
{
  PRECISION u[ 3], q[ 3];
  PRECISION omega, f;
  PRECISION K, h, g;

#define matrix( i, j) MATRIX[ i * 4 + j]
#define eigen_vectors( i, j) EIGEN_VECTORS[ i * 4 + j]

  // Initialize EIGEN_VECTORS to the identity matrix
  for( size_t i = 0; i < 3; ++i)
  {
    eigen_vectors( i, i) = 1.0;
    for( size_t j = 0; j < i; ++j)
    {
      eigen_vectors( i, j) = eigen_vectors( j, i) = 0.0;
    }
  }

  // Bring first row and column to the desired form
  h = pown( matrix( 0, 1), 2) + pown( matrix( 0, 2), 2);
  if( matrix( 0, 1) > 0.0)
  {
    g = -sqrt( h);
  }
  else
  {
    g = sqrt( h);
  }
  OFF_DIAGONAL[ 0] = g;
  f     = g * matrix( 0, 1);
  u[ 1] = matrix( 0, 1) - g;
  u[ 2] = matrix( 0, 2);

  omega = h - f;
  if( omega > 0.0)
  {
    omega = 1.0 / omega;
    K     = 0.0;
    for( size_t i = 1; i < 3; ++i)
    {
      f      = matrix( 1, i) * u[ 1] + matrix( i, 2) * u[2];
      q[ i]  = omega * f;                  // p
      K     += u[ i] * f;                  // u* A u
    }
    K *= 0.5 * pown( omega, 2);

    for( size_t i = 1; i < 3; ++i)
    {
      q[ i] = q[ i] - K * u[ i];
    }

    EIGENVALUES[ 0] = matrix( 0, 0);
    EIGENVALUES[ 1] = matrix( 1, 1) - 2 * q[ 1] * u[ 1];
    EIGENVALUES[ 2] = matrix( 2, 2) - 2 * q[ 2] * u[ 2];

    // Store inverse Householder transformation in Q
    for( size_t j = 1; j < 3; ++j)
    {
      f = omega * u[ j];
      for( size_t i = 1; i < 3; ++i)
      {
        eigen_vectors( i, j) -= f * u[ i];
      }
    }

    // Calculate updated A[1][2] and store it in OFF_DIAGONAL[1]
    OFF_DIAGONAL[ 1] = matrix( 1, 2) - q[ 1] * u[ 2] - u[ 1] * q[ 2];
  }
  else
  {
    for( size_t i = 0; i < 3; ++i)
    {
      EIGENVALUES[ i] = matrix( i, i);
    }
    OFF_DIAGONAL[ 1] = matrix( 1, 2);
  }
#undef matrix
#undef eigen_vectors
}

//! @brief Eigenvalues and Eigenvectors of a 3x3 matrix using ql method
//! @param MATRIX the matrix of size localsize( 0) x localsize( 1)
//! @param EIGENVECTORS matrix to write eigenvectors to
//! @param EIGENVALUES the eigenvalue vector
void Matrix3x3SymmetricEigenvectorsQLAux
(
  const PRECISION *MATRIX,
        PRECISION *EIGEN_VECTORS,
        PRECISION *EIGEN_VALUES
)
{
  PRECISION off_diagonal[ 3];   // The third element is used only as temporary workspace
  PRECISION g, r, p, f, b, s, c, t; // Intermediate storage

  // Transform A to real tridiagonal form by the Householder method
  // copy to shared memory
#define eigen_vectors( i, j) EIGEN_VECTORS[ i * 4 + j]

  Matrix3x3SymmetricTridiagonalizeHouseholderAux( MATRIX, EIGEN_VECTORS, EIGEN_VALUES, off_diagonal);

  // Calculate eigensystem of the remaining real symmetric tridiagonal matrix
  // with the QL method
  //
  // Loop over all off-diagonal elements
  for( int l = 0; l < 3 - 1; ++l)
  {
    uint num_iter = 0;
    while ( true)
    {
      // Check for convergence and exit iteration loop if off-diagonal
      // element e(l) is zero
      int m;
      for( m = l; m <= 3 - 2; ++m)
      {
        g = fabs( EIGEN_VALUES[ m]) + fabs( EIGEN_VALUES[ m + 1]);
        if( fabs( off_diagonal[ m]) + g == g)
        {
          break;
        }
      }
      if( m == l)
      {
        break;
      }

      // should converge within 30 iterations
      if( num_iter++ >= 30)
      {
        return;
      }

      // Calculate g = d_m - k
      g = ( EIGEN_VALUES[ l + 1] - EIGEN_VALUES[ l]) / ( off_diagonal[ l] + off_diagonal[ l]);
      r = sqrt( pown( g, 2) + 1.0);
      if( g > 0.0)
      {
        g = EIGEN_VALUES[ m] - EIGEN_VALUES[ l] + off_diagonal[ l] / ( g + r);
      }
      else
      {
        g = EIGEN_VALUES[ m] - EIGEN_VALUES[ l] + off_diagonal[ l] / ( g - r);
      }

      s = c = 1.0;
      p = 0.0;
      for( int i = m - 1; i >= l; --i)
      {
        f = s * off_diagonal[ i];
        b = c * off_diagonal[ i];
        if( fabs( f) > fabs( g))
        {
          c      = g / f;
          r      = sqrt( pown( c, 2) + 1.0);
          off_diagonal[ i + 1] = f * r;
          c     *= (s = 1.0/r);
        }
        else
        {
          s      = f / g;
          r      = sqrt( pown( s, 2) + 1.0);
          off_diagonal[ i + 1] = g * r;
          s     *= ( c = 1.0 / r);
        }

        g = EIGEN_VALUES[ i + 1] - p;
        r = ( EIGEN_VALUES[ i] - g) * s + 2 * c * b;
        p = s * r;
        EIGEN_VALUES[ i + 1] = g + p;
        g = c * r - b;

        // Form eigenvectors
        for( size_t k = 0; k < 3; ++k)
        {
          t = eigen_vectors( k, i + 1);
          eigen_vectors( k, i + 1) = s * eigen_vectors( k, i) + c * t;
          eigen_vectors( k, i    ) = c * eigen_vectors( k, i) - s * t;
        }
      }
      EIGEN_VALUES[ l] -= p;
      off_diagonal[ l] = g;
      off_diagonal[ m] = 0.0;
    }
  }
#undef eigen_vectors
}

//! @brief Eigenvalues and Eigenvectors of a 3x3 matrix
//! @param MATRIX the matrix of size localsize( 0) x localsize( 1)
//! @param EIGENVECTORS matrix to write eigenvectors to
//! @param EIGENVALUES the eigenvalue vector
void Matrix3x3SymmetricEigenvectorsAux
(
  const PRECISION *MATRIX,
        PRECISION *EIGEN_VECTORS,
        PRECISION *EIGEN_VALUES
)
{
  // calculate eigenvalues
  Matrix3x3SymmetricEigenvaluesAux( MATRIX, 0, EIGEN_VALUES);

#define matrix( i, j) MATRIX[ i * 4 + j]
#define eigen_vectors( i, j) EIGEN_VECTORS[ i * 4 + j]

  const PRECISION n0 = pown( matrix( 0, 0), 2) + pown( matrix( 0, 1), 2) + pown( matrix( 0, 2), 2);
  const PRECISION n1 = pown( matrix( 0, 1), 2) + pown( matrix( 1, 1), 2) + pown( matrix( 1, 2), 2);

  PRECISION u, norm;
  PRECISION t = fabs( EIGEN_VALUES[ 0]);
  if( ( u = fabs( EIGEN_VALUES[ 1])) > t)
  {
    t = u;
  }
  if( ( u = fabs( EIGEN_VALUES[ 2])) > t)
  {
    t = u;
  }
  if( t < 1.0)
  {
    u = t;
  }
  else
  {
    u = pown( t, 2);
  }
  PRECISION error = 256 * FLT_EPSILON * ( n0 + u) * ( n1 + u);

  eigen_vectors( 0, 1) = matrix( 0, 1) * matrix( 1, 2) - matrix( 0,2) * matrix( 1, 1);
  eigen_vectors( 1, 1) = matrix( 0, 2) * matrix( 0, 1) - matrix( 1,2) * matrix( 0, 0);
  eigen_vectors( 2, 1) = pown( matrix( 0, 1), 2);

  // Calculate first eigenvector by the formula
  //   v[0] = (A - EIGEN_VALUES( 0)).e1 x (A - EIGEN_VALUES( 0)).e2
  eigen_vectors( 0, 0) = eigen_vectors( 0, 1) + matrix( 0, 2) * EIGEN_VALUES[ 0];
  eigen_vectors( 1, 0) = eigen_vectors( 1, 1) + matrix( 1, 2) * EIGEN_VALUES[ 0];
  eigen_vectors( 2, 0) = ( matrix( 0, 0) - EIGEN_VALUES[ 0]) * ( matrix( 1, 1) - EIGEN_VALUES[ 0]) - eigen_vectors( 2, 1);
  norm = pown( eigen_vectors( 0, 0), 2) + pown( eigen_vectors( 1, 0), 2) + pown( eigen_vectors( 2, 0), 2);

  // If vectors are nearly linearly dependent, or if there might have
  // been large cancellations in the calculation of A[i][i] - EIGEN_VALUES( 0), fall
  // back to QL algorithm
  // Note that this simultaneously ensures that multiple eigenvalues do
  // not cause problems: If EIGEN_VALUES( 0) = EIGEN_VALUES( 1), then A - EIGEN_VALUES( 0) * I has rank 1,
  // i.e. all columns of A - EIGEN_VALUES( 0) * I are linearly dependent.
  if( norm <= error)
  {
    Matrix3x3SymmetricEigenvectorsQLAux( MATRIX, EIGEN_VECTORS, EIGEN_VALUES);
    return;
  }

  norm = pown( 1.0 / norm, 2);
  for( size_t j = 0; j < 3; ++j)
  {
    eigen_vectors( j, 0) *= norm;
  }

  // calculate second eigenvector by the formula
  //   v[1] = (A - EIGEN_VALUES( 1)).e1 x (A - EIGEN_VALUES( 1)).e2
  eigen_vectors( 0, 1) = eigen_vectors( 0, 1) + matrix( 0, 2) * EIGEN_VALUES[ 1];
  eigen_vectors( 1, 1) = eigen_vectors( 1, 1) + matrix( 1, 2) * EIGEN_VALUES[ 1];
  eigen_vectors( 2, 1) = matrix( 0, 0) * matrix( 1, 1) - eigen_vectors( 2, 1);
  norm  = pown( eigen_vectors( 0, 1), 2) + pown( eigen_vectors( 1, 1), 2) + pown( eigen_vectors( 2, 1), 2);

  error = n0 * n1;
  if( norm <= error)
  {
    Matrix3x3SymmetricEigenvectorsQLAux( MATRIX, EIGEN_VECTORS, EIGEN_VALUES);
    return;
  }

  norm = pown( 1.0 / norm, 2);
  for( size_t j = 0; j < 3; ++j)
  {
    eigen_vectors( j, 1) *= norm;
  }

  // Calculate third eigenvector according to
  //   v[2] = v[0] x v[1]
  eigen_vectors( 0, 2) = eigen_vectors( 1, 0) * eigen_vectors( 2, 1) - eigen_vectors( 2,0) * eigen_vectors( 1, 1);
  eigen_vectors( 1, 2) = eigen_vectors( 2, 0) * eigen_vectors( 0, 1) - eigen_vectors( 0,0) * eigen_vectors( 2, 1);
  eigen_vectors( 2, 2) = eigen_vectors( 0, 0) * eigen_vectors( 1, 1) - eigen_vectors( 1,0) * eigen_vectors( 0, 1);

#undef shared_matrix
#undef eigen_vectors
}

//! @brief Eigenvalues and Eigenvectors of a 3x3 matrix
//! @param MATRIX the matrix of size localsize( 0) x localsize( 1)
//! @param EIGENVECTORS matrix to write eigenvectors to
//! @param EIGENVALUES the eigenvalue vector
__kernel void Matrix3x3SymmetricEigenvectors
(
  const __global PRECISION *MATRIX,
        __global PRECISION *EIGEN_VECTORS,
        __global PRECISION *EIGEN_VALUES
)
{
  PRECISION matrix_local[ 16];
  PRECISION eigen_vectors_local[ 16];
  PRECISION eigen_values_local[ 4];

  // copy input matrix
  for( uint row = 0; row < 4; ++row)
  {
    eigen_values_local[ row] = 0.0;
    for( uint col = 0; col < 4; ++col)
    {
      eigen_vectors_local[ row * 4 + col] = 0.0;
      matrix_local[ row * 4 + col] = MATRIX[ row * 16 + col];
    }
  }

  // calculate with auxilary functions
  Matrix3x3SymmetricEigenvectorsAux( matrix_local, eigen_vectors_local, eigen_values_local);

  // copy results
  for( uint row = 0; row < 4; ++row)
  {
    EIGEN_VALUES[ row] = eigen_values_local[ row];
    for( uint col = 0; col < 4; ++col)
    {
      EIGEN_VECTORS[ row * 16 + col] = eigen_vectors_local[ row * 4 + col];
    }
  }
}

//! @brief sort rows by eigenvalues
//! @param MATRIX the matrix of size localsize(0) x localsize(0)
//! @param EIGEN_VALUES the eigenvalue vector
__kernel void Matrix3x3SortRowsAndVector
(
  __global PRECISION *MATRIX,
  __global PRECISION *EIGEN_VALUES,
  __local  PRECISION *SHARED_M,
  __local  PRECISION *SHARED_V
)
{
  const uint block_size = get_local_size( 0);
  const uint thread_id =  get_local_id( 0);

  // copy elements row-wise (each thread one row
  SHARED_V[ thread_id] = EIGEN_VALUES[ thread_id];
  for( uint col = 0; col < block_size; ++col)
  {
    SHARED_M[ thread_id * block_size + col] = MATRIX[ thread_id * block_size + col];
  }

  barrier( CLK_LOCAL_MEM_FENCE);
  PRECISION tmp;

  if( SHARED_V[ 0] < SHARED_V[ 1])
  {
    barrier( CLK_LOCAL_MEM_FENCE);
    if( thread_id == 0)
    {
      tmp = SHARED_V[ 0];
      SHARED_V[ 0] = SHARED_V[ 1];
      SHARED_V[ 1] = tmp;
    }
    tmp = SHARED_M[ 0 * block_size + thread_id];
    SHARED_M[ 0 * block_size + thread_id] = SHARED_M[ 1 * block_size + thread_id];
    SHARED_M[ 1 * block_size + thread_id] = tmp;
  }

  barrier( CLK_LOCAL_MEM_FENCE);

  if( SHARED_V[ 1] < SHARED_V[ 2])
  {
    barrier( CLK_LOCAL_MEM_FENCE);
    if( thread_id == 0)
    {
      tmp = SHARED_V[ 1];
      SHARED_V[ 1] = SHARED_V[ 2];
      SHARED_V[ 2] = tmp;
    }
    tmp = SHARED_M[ 1 * block_size + thread_id];
    SHARED_M[ 1 * block_size + thread_id] = SHARED_M[ 2 * block_size + thread_id];
    SHARED_M[ 2 * block_size + thread_id] = tmp;
  }

  barrier( CLK_LOCAL_MEM_FENCE);

  if( SHARED_V[ 0] < SHARED_V[ 1])
  {
    barrier( CLK_LOCAL_MEM_FENCE);
    if( thread_id == 0)
    {
      tmp = SHARED_V[ 0];
      SHARED_V[ 0] = SHARED_V[ 1];
      SHARED_V[ 1] = tmp;
    }
    tmp = SHARED_M[ 0 * block_size + thread_id];
    SHARED_M[ 0 * block_size + thread_id] = SHARED_M[ 1 * block_size + thread_id];
    SHARED_M[ 1 * block_size + thread_id] = tmp;
  }

  barrier( CLK_LOCAL_MEM_FENCE);

  EIGEN_VALUES[ thread_id] = SHARED_V[ thread_id];
  for( uint row = 0; row < block_size; ++row)
  {
    MATRIX[ row * block_size + thread_id] = SHARED_M[ row * block_size + thread_id];
  }
}

//! @brief replace the Row with the cross product of the other two rows
//! @param MATRIX 3x3 matrix
//! @param ROW the row to orthogonalize by the cross product of the other two
__kernel void Matrix3x3OrthogonalizeRow
(
        __global PRECISION *MATRIX,
  const          uint       ROW
)
{
  const size_t local_index_x = get_local_id( 0);
  const size_t block_size = get_local_size( 0);
  if( local_index_x >= 3)
  {
    return;
  }
  MATRIX[ ROW * block_size + local_index_x] =
      MATRIX[ ( ( ROW + 1) % 3) * block_size + ( local_index_x + 1) % 3] * MATRIX[ ( ( ROW + 2) % 3) * block_size + ( local_index_x + 2) % 3] -
      MATRIX[ ( ( ROW + 1) % 3) * block_size + ( local_index_x + 2) % 3] * MATRIX[ ( ( ROW + 2) % 3) * block_size + ( local_index_x + 1) % 3];
}

//! @brief normalize rows by square of corresponding elements in vector
//! @param MATRIX the matrix of interest
//! @param VECTOR the vector of normalization values
__kernel void Matrix3x3NormalizeRows
(
        __global PRECISION *MATRIX,
  const __global PRECISION *VECTOR
)
{
  const size_t local_index_x = get_local_id( 0);
  const size_t local_index_y = get_local_id( 1);
  const size_t block_size    = get_local_size( 0);
  const PRECISION value = 1.0 / sqrt( VECTOR[ local_index_y]);
  if( isnan( value) || isinf( value))
  {
    return;
  }
  MATRIX[ local_index_y * block_size + local_index_x] *= value;
}

//! @brief multiply matrix with its transposed
//! @param MATRIX
//! @param SHARED_MATRIX
__kernel void Matrix3x3MultiplyWithTransposed
(
  __global PRECISION* MATRIX,
  __local  PRECISION* SHARED_MATRIX
)
{
  const size_t size_x = get_local_size( 0);

  const size_t index_x = get_local_id( 0);
  const size_t index_y = get_local_id( 1);

  SHARED_MATRIX[ index_y * size_x + index_x] = MATRIX[ index_y * size_x + index_x];
  barrier( CLK_LOCAL_MEM_FENCE);

  PRECISION accumulator = 0;

  for( size_t i = 0; i < size_x; ++i)
  {
    accumulator += SHARED_MATRIX[ index_y * size_x + i] * SHARED_MATRIX[ index_x * size_x + i];
  }

  MATRIX[ index_y * size_x + index_x] = accumulator;
}

//! @brief apply translation to transformation matrix
//! @param TRANSFORMATION_MATRIX matrix with rotation and transformation
//! @param TRANSLATION_VECTOR vector of translation
//! @param SCALING multiply TRANSLATION by the scaling
__kernel void TranslationOnTransformation
(
        __global PRECISION *TRANSFORMATION_MATRIX,
  const __global PRECISION *TRANSLATION_VECTOR,
  const          PRECISION  SCALING
)
{
  const size_t local_index_x = get_local_id( 0);
  const size_t block_size = get_local_size( 0);

  // add translation to third row
  TRANSFORMATION_MATRIX[ 3 * block_size + local_index_x] += SCALING * TRANSLATION_VECTOR[ local_index_x];
}

//! @brief apply rotation to transformation matrix
//! @param TRANSFORMATION_MATRIX matrix with rotation and transformation
//! @param ROTATION_MATRIX rotation matrix
__kernel void RotationOnTransformation
(
        __global PRECISION *TRANSFORMATION_MATRIX,
  const __global PRECISION *ROTATION_MATRIX,
        __local  PRECISION *SHARED_TRANS,
        __local  PRECISION *SHARED_ROT
)
{
  const size_t local_index_x = get_local_id( 1);
  const size_t local_index_y = get_local_id( 0);
  const size_t block_size = get_local_size( 0);

  // copy to shared memory
  SHARED_TRANS[ local_index_y * block_size + local_index_x] = TRANSFORMATION_MATRIX[ local_index_y * block_size + local_index_x];
  SHARED_ROT[ local_index_y * block_size + local_index_x] = ROTATION_MATRIX[ local_index_y * block_size + local_index_x];

  // set scaling factor for rotation matrix, so that translation and scaling does not get lost
  if( local_index_x == 3 && local_index_y == 3)
  {
    SHARED_ROT[ local_index_y * block_size + local_index_x] = 1.0;
  }

  // sync shared memory copy
  barrier( CLK_LOCAL_MEM_FENCE);

  // multiplication
  PRECISION accumulator = 0.0;

  for( size_t i = 0; i < 4; ++i)
  {
    accumulator += SHARED_TRANS[ local_index_y * block_size + i] * SHARED_ROT[ i * block_size + local_index_x];
  }

  // set matrix
  TRANSFORMATION_MATRIX[ local_index_y * block_size + local_index_x] = accumulator;
}

__kernel void TransformationsFromCovarianceMatrices
(
         __global PRECISION *COVARIANCE_MATRICES,
   const __global PRECISION *CENTER_COORDINATES,
   const __global PRECISION *CENTER_REFERENCE_COORDINATES,
   const               uint  NUMBER_COLS,
   const               uint  NUMBER_ROWS_MATRIX,
   const               uint  NUMBER_MATRICES
)
{
  const size_t mat_num = get_global_id( 0);

  if( mat_num >= NUMBER_MATRICES)
  {
    return;
  }

  PRECISION moment_matrix[ 16];
  PRECISION rotate_matrix[ 16];
  PRECISION eigenvectors[  16];
  PRECISION tmp_rot[       16];
  for( size_t i = 0; i < 16; ++i)
  {
    moment_matrix[ i] = 0.0;
    eigenvectors[  i] = 0.0;
    rotate_matrix[ i] = 0.0;
    tmp_rot[       i] = 0.0;
  }

  // copy from global memory
  for( size_t row = 0; row < 3; ++row)
  {
    for( size_t col = 0; col < 3; ++col)
    {
      moment_matrix[ row * 4 + col] = COVARIANCE_MATRICES[ mat_num * NUMBER_ROWS_MATRIX * NUMBER_COLS + row * NUMBER_COLS + col];
    }
  }

  // multiply with its transposed to symmetrize
  for( size_t row = 0; row < 3; ++row)
  {
    for( size_t col = 0; col < 3; ++col)
    {
      for( size_t i = 0; i < 3; ++i)
      {
        rotate_matrix[ row * 4 + col] += moment_matrix[ row * 4 + i] * moment_matrix[ col * 4 + i];
      }
    }
  }
  // determine eigenvalues and eigenvectors
  PRECISION eigenvalues[ 3];
  Matrix3x3SymmetricEigenvectorsAux( rotate_matrix, eigenvectors, eigenvalues);

  // sort cols by eigenvalues (on cpu it is transposed, sort rows and transposed back)
  PRECISION tmp;

  if( eigenvalues[ 0] < eigenvalues[ 1])
  {
    tmp = eigenvalues[ 0];
    eigenvalues[ 0] = eigenvalues[ 1];
    eigenvalues[ 1] = tmp;

    for( size_t i = 0; i < 3; ++ i)
    {
      tmp = eigenvectors[ i * 4 + 0];
      eigenvectors[ i * 4 + 0] = eigenvectors[ i * 4 + 1];
      eigenvectors[ i * 4 + 1] = tmp;
    }
  }

  if( eigenvalues[ 1] < eigenvalues[ 2])
  {
    tmp = eigenvalues[ 1];
    eigenvalues[ 1] = eigenvalues[ 2];
    eigenvalues[ 2] = tmp;

    for( size_t i = 0; i < 3; ++ i)
    {
      tmp = eigenvectors[ i * 4 + 1];
      eigenvectors[ i * 4 + 1] = eigenvectors[ i * 4 + 2];
      eigenvectors[ i * 4 + 2] = tmp;
    }
  }

  if( eigenvalues[ 0] < eigenvalues[ 1])
  {
    tmp = eigenvalues[ 0];
    eigenvalues[ 0] = eigenvalues[ 1];
    eigenvalues[ 1] = tmp;

    for( size_t i = 0; i < 3; ++ i)
    {
      tmp = eigenvectors[ i * 4 + 0];
      eigenvectors[ i * 4 + 0] = eigenvectors[ i * 4 + 1];
      eigenvectors[ i * 4 + 1] = tmp;
    }
  }

  // orthogonalize last col
  eigenvectors[ 0 * 4 + 2] = eigenvectors[ 1 * 4 + 0] * eigenvectors[ 2 * 4 + 1] - eigenvectors[ 2 * 4 + 0] * eigenvectors[ 1 * 4 + 1];
  eigenvectors[ 1 * 4 + 2] = eigenvectors[ 2 * 4 + 0] * eigenvectors[ 0 * 4 + 1] - eigenvectors[ 0 * 4 + 0] * eigenvectors[ 2 * 4 + 1];
  eigenvectors[ 2 * 4 + 2] = eigenvectors[ 0 * 4 + 0] * eigenvectors[ 1 * 4 + 1] - eigenvectors[ 1 * 4 + 0] * eigenvectors[ 0 * 4 + 1];

  // multiply by moment but remember that rotate isn't transposed as it should be
  for( size_t i = 0; i < 3; ++i)
  {
    for( size_t j = 0; j < 3; ++j)
    {
      rotate_matrix[ i * 4 + j] = 0.0;
      for( size_t k = 0; k < 4; ++k)
      {
        rotate_matrix[ i * 4 + j] += eigenvectors[ k * 4 + i] * moment_matrix[ k * 4 + j];
      }
    }
  }

  // normalize columns
  for( size_t row = 0; row < 3; ++row)
  {
    const PRECISION value = 1.0 / sqrt( eigenvalues[ row]);
//    if( isnan( value) || isinf( value))
//    {
//      continue;
//    }
    for( size_t col = 0; col < 3; ++col)
    {
      rotate_matrix[ row * 4 + col] *= value;
    }
  }

  // orthogonalize last row again
  rotate_matrix[ 2 * 4 + 0] = rotate_matrix[ 0 * 4 + 1] * rotate_matrix[ 1 * 4 + 2] - rotate_matrix[ 0 * 4 + 2] * rotate_matrix[ 1 * 4 + 1];
  rotate_matrix[ 2 * 4 + 1] = rotate_matrix[ 0 * 4 + 2] * rotate_matrix[ 1 * 4 + 0] - rotate_matrix[ 0 * 4 + 0] * rotate_matrix[ 1 * 4 + 2];
  rotate_matrix[ 2 * 4 + 2] = rotate_matrix[ 0 * 4 + 0] * rotate_matrix[ 1 * 4 + 1] - rotate_matrix[ 0 * 4 + 1] * rotate_matrix[ 1 * 4 + 0];

  // multiply by eigenvectors - both eigenvectors and rotate should have been transposed and weren't
  for( size_t i = 0; i < 3; ++i)
  {
    for( size_t j = 0; j < 3; ++j)
    {
      for( size_t k = 0; k < 3; ++k)
      {
        tmp_rot[ i * 4 + j] += rotate_matrix[ k * 4 + i] * eigenvectors[ j * 4 + k];
      }
    }
  }
  tmp_rot[ 3 * 4 + 3] = 1.0;

  for( size_t i = 0; i < 16; ++i)
  {
    rotate_matrix[ i] = 0.0;
    moment_matrix[ i] = 0.0;
  }

  // set diagonal of transformation matrix to 1
  for( size_t i = 0; i < 4; ++i)
  {
    moment_matrix[ i * 4 + i] = 1.0;
  }

  PRECISION ref_coord_center[ 3], coord_center[ 3];
  for( size_t i = 0; i < 3; ++i)
  {
    ref_coord_center[ i] = CENTER_REFERENCE_COORDINATES[ mat_num * NUMBER_COLS + i];
    coord_center[     i] = CENTER_COORDINATES[           mat_num * NUMBER_COLS + i];
  }

  // apply translation on transformation matrix which was previously a separate kernel "TranslationOnTransformation"
  // add translation to third row
  for( size_t i = 0; i < 3; ++i)
  {
    moment_matrix[ 3 * 4 + i] -= coord_center[ i];
  }

  // apply rotation on transformation matrix which was also a kernel "RotationOnTransformation"
  for( size_t i = 0; i < 4; ++i)
  {
    for( size_t j = 0; j < 4; ++j)
    {
      for( size_t k = 0; k < 4; ++k)
      {
        rotate_matrix[ i * 4 + j] += moment_matrix[ i * 4 + k] * tmp_rot[ k * 4 + j];
      }
    }
  }

  // apply translation on transformation matrix which was previously a separate kernel "TranslationOnTransformation"
  // add translation to third row
  for( size_t i = 0; i < 3; ++i)
  {
    rotate_matrix[ 3 * 4 + i] += ref_coord_center[ i];
  }

  // return transformation matrix
  for( size_t row = 0; row < NUMBER_ROWS_MATRIX; ++row)
  {
    for( size_t col = 0; col < 4; ++col)
    {
      COVARIANCE_MATRICES[ mat_num * NUMBER_ROWS_MATRIX * NUMBER_COLS + row * NUMBER_COLS + col] = rotate_matrix[ row * 4 + col];
    }
  }
}
