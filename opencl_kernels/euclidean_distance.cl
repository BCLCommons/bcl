// this is only included, so that eclipse c++ editor can be used and parses information correctly
// this will not be interpreted, when the program is compiled with an OpenCL compiler
#ifndef __OPENCL_VERSION__
 #include "kernel.h"
#endif

//! @brief kernel for calculating pairwise distances between rows of two matrices yielding distance matrix
//! @param OUTPUT buffer to place resulting distance matrix
//! @param INPUT_A, INPUT_B the two input matrices
//! @param As, Bs the local memory allocations
//! @param ROWS_A, COLS_A, ROWS_B, COLS_B the dimensionality
__kernel void EuclideanDistance
(
  __global PRECISION* OUTPUT,
  __const __global PRECISION* INPUT_A,
  __const __global PRECISION* INPUT_B,
  __local PRECISION* As,
  __local PRECISION* Bs,
  __const int ROWS_A,
  __const int COLS_A,
  __const int ROWS_B,
  __const int COLS_B
)
{
  int bx = get_group_id( 0), by = get_group_id( 1);
  int tx = get_local_id( 0), ty = get_local_id( 1);
  int block_size = get_local_size( 0);
  int A_Begin = by * block_size * COLS_A;
  int B_Begin = bx * block_size * COLS_B;
  int A_End = A_Begin + COLS_A - 1, A_idx, B_idx, k, o;
  PRECISION tmp, s = 0.0;

  for( A_idx = A_Begin, B_idx = B_Begin; A_idx <= A_End; A_idx += block_size, B_idx += block_size)
  {
    As[ ty * block_size + tx] = INPUT_A[ A_idx + ty * COLS_A + tx];
    Bs[ tx * block_size + ty] = INPUT_B[ B_idx + ty * COLS_B + tx];
    barrier( CLK_LOCAL_MEM_FENCE);
//#pragma unroll
    for( k = 0; k < block_size; k++)
    {
      tmp = As[ ty * block_size + k] - Bs[ k * block_size + tx];
      s += tmp * tmp;
    }
    barrier( CLK_LOCAL_MEM_FENCE);
  }

  o = by * block_size * ROWS_A + ty * ROWS_B + bx * block_size + tx;
  OUTPUT[ o] = sqrt( s);
}

//! @brief kernel for calculating pairwise distances between rows of two matrices yielding distance matrix
//! @param OUTPUT buffer to place resulting distance matrix
//! @param INPUT_A, INPUT_B the two input matrices
//! @param As, Bs the local memory allocations
//! @param ROWS_A, COLS_A, ROWS_B, COLS_B the dimensionality
__kernel void EuclideanDistanceFast
(
  __global float* OUTPUT,
  __const __global float* INPUT_A,
  __const __global float* INPUT_B,
  __local float* As,
  __local float* Bs,
  __const int ROWS_A,
  __const int COLS_A,
  __const int ROWS_B,
  __const int COLS_B
)
{
  int bx = get_group_id( 0), by = get_group_id( 1);
  int tx = get_local_id( 0), ty = get_local_id( 1);
  int block_size = get_local_size( 0);
  int A_Begin = by * block_size * COLS_A;
  int B_Begin = bx * block_size * COLS_B;
  int A_End = A_Begin + COLS_A - 1, A_idx, B_idx, o;
  float4 s = ( float4)( 0.0, 0.0, 0.0, 0.0);

  for( A_idx = A_Begin, B_idx = B_Begin; A_idx <= A_End; A_idx += block_size, B_idx += block_size)
  {
    As[ ty * block_size + tx] = INPUT_A[ A_idx + ty * COLS_A + tx];
    Bs[ tx * block_size + ty] = INPUT_B[ B_idx + ty * COLS_B + tx];
    barrier( CLK_LOCAL_MEM_FENCE);

    float4 tmp_a0 = ( float4)( As[ ty * block_size + 0], As[ ty * block_size + 1], As[ ty * block_size + 2], As[ ty * block_size + 3]);
    float4 tmp_b0 = ( float4)( Bs[ 0 * block_size + tx], Bs[ 1 * block_size + tx], Bs[ 2 * block_size + tx], Bs[ 3 * block_size + tx]);
    float4 tmp_c0 =  tmp_a0 - tmp_b0;
    s = mad( tmp_c0, tmp_c0, s);

    float4 tmp_a1 = ( float4)( As[ ty * block_size + 4], As[ ty * block_size + 5], As[ ty * block_size + 6], As[ ty * block_size + 7]);
    float4 tmp_b1 = ( float4)( Bs[ 4 * block_size + tx], Bs[ 5 * block_size + tx], Bs[ 6 * block_size + tx], Bs[ 7 * block_size + tx]);
    float4 tmp_c1 =  tmp_a1 - tmp_b1;
    s = mad( tmp_c1, tmp_c1, s);

    float4 tmp_a2 = ( float4)( As[ ty * block_size + 8], As[ ty * block_size + 9], As[ ty * block_size + 10], As[ ty * block_size + 11]);
    float4 tmp_b2 = ( float4)( Bs[ 8 * block_size + tx], Bs[ 9 * block_size + tx], Bs[ 10 * block_size + tx], Bs[ 11 * block_size + tx]);
    float4 tmp_c2 =  tmp_a2 - tmp_b2;
    s = mad( tmp_c2, tmp_c2, s);

    float4 tmp_a3 = ( float4)( As[ ty * block_size + 12], As[ ty * block_size + 13], As[ ty * block_size + 14], As[ ty * block_size + 15]);
    float4 tmp_b3 = ( float4)( Bs[ 12 * block_size + tx], Bs[ 13 * block_size + tx], Bs[ 14 * block_size + tx], Bs[ 15 * block_size + tx]);
    float4 tmp_c3 =  tmp_a3 - tmp_b3;
    s = mad( tmp_c3, tmp_c3, s);

    barrier( CLK_LOCAL_MEM_FENCE);
  }
  float4 in = ( float4)( 1, 1, 1, 1);
  o = by * block_size * ROWS_A + ty * ROWS_B + bx * block_size + tx;
  OUTPUT[ o] = native_sqrt( dot( in, s));
}
