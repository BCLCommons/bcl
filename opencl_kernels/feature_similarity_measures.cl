// this is only included, so that eclipse c++ editor can be used and parses information correctly
// this will not be interpreted, when the program is compiled with an OpenCL compiler
#ifndef __OPENCL_VERSION__
 #include "kernel.h"
#endif

#define BLOCK_SIZE 16

//! @brief kernel for calculating pairwise euclidean similarity coefficients
//! @param OUTPUT buffer to place resulting distance matrix
//! @param INPUT_A, INPUT_A the two input matrices
//! @param As, Bs the local memory allocations
//! @param ROWS, COLS, ROWS, COLS the dimensionality
__kernel void Euclidean
(
  __global PRECISION* OUTPUT,
  __const __global PRECISION* INPUT_A,
  __local PRECISION* As,
  __local PRECISION* Bs,
  __const int ROWS,
  __const int COLS
)
{
  int bx = get_group_id( 0), by = get_group_id( 1);
  int tx = get_local_id( 0), ty = get_local_id( 1);
  int A_Begin = by * BLOCK_SIZE * COLS;
  int B_Begin = bx * BLOCK_SIZE * COLS;
  int A_End = A_Begin + COLS - 1, A_idx, B_idx, k, o;
  PRECISION tmp, s = 0.0;

  for( A_idx = A_Begin, B_idx = B_Begin; A_idx <= A_End; A_idx += BLOCK_SIZE, B_idx += BLOCK_SIZE)
  {
    As[ ty * BLOCK_SIZE + tx] = INPUT_A[ A_idx + ty * COLS + tx];
    Bs[ tx * BLOCK_SIZE + ty] = INPUT_A[ B_idx + ty * COLS + tx];
    barrier( CLK_LOCAL_MEM_FENCE);
#pragma unroll
    for( k = 0; k < BLOCK_SIZE; k++)
    {
      tmp = As[ ty * BLOCK_SIZE + k] - Bs[ k * BLOCK_SIZE + tx];
      s += tmp * tmp;
    }
    barrier( CLK_LOCAL_MEM_FENCE);
  }

  o = by * BLOCK_SIZE * ROWS + ty * ROWS + bx * BLOCK_SIZE + tx;
  OUTPUT[ o] = native_recip( 1 + native_sqrt( s));
}


//! @brief kernel for calculating pairwise tanimoto similarity coefficients
//! @param OUTPUT buffer to place resulting distance matrix
//! @param INPUT_A, INPUT_A the two input matrices
//! @param As, Bs the local memory allocations
//! @param ROWS, COLS, ROWS, COLS the dimensionality
__kernel void Tanimoto
(
  __global PRECISION* OUTPUT,
  __const __global PRECISION* INPUT_A,
  __local PRECISION* As,
  __local PRECISION* Bs,
  __const int ROWS,
  __const int COLS
)
{
  int bx = get_group_id( 0), by = get_group_id( 1);
  int tx = get_local_id( 0), ty = get_local_id( 1);
  int A_Begin = by * BLOCK_SIZE * COLS;
  int B_Begin = bx * BLOCK_SIZE * COLS;
  int A_End = A_Begin + COLS - 1, A_idx, B_idx, k, o;
  PRECISION tmp, a_b_sum = 0.0, a2_sum = 0.0, b2_sum = 0.0, a = 0.0, b = 0.0;

  for( A_idx = A_Begin, B_idx = B_Begin; A_idx <= A_End; A_idx += BLOCK_SIZE, B_idx += BLOCK_SIZE)
  {
    As[ ty * BLOCK_SIZE + tx] = INPUT_A[ A_idx + ty * COLS + tx];
    Bs[ tx * BLOCK_SIZE + ty] = INPUT_A[ B_idx + ty * COLS + tx];
    barrier( CLK_LOCAL_MEM_FENCE);
#pragma unroll
    for( k = 0; k < BLOCK_SIZE; k++)
    {
      a = As[ ty * BLOCK_SIZE + k];
      b = Bs[ k * BLOCK_SIZE + tx];
      a_b_sum += a * b;
      a2_sum  += a * a;
      b2_sum  += b * b;
    }
    barrier( CLK_LOCAL_MEM_FENCE);
  }

  o = by * BLOCK_SIZE * ROWS + ty * ROWS + bx * BLOCK_SIZE + tx;
  OUTPUT[ o] = a_b_sum * native_recip( a2_sum + b2_sum - a_b_sum);
}

//! @brief kernel for calculating pairwise Dice similarity coefficients
//! @param OUTPUT buffer to place resulting distance matrix
//! @param INPUT_A, INPUT_A the two input matrices
//! @param As, Bs the local memory allocations
//! @param ROWS, COLS, ROWS, COLS the dimensionality
__kernel void Dice
(
  __global PRECISION* OUTPUT,
  __const __global PRECISION* INPUT_A,
  __local PRECISION* As,
  __local PRECISION* Bs,
  __const int ROWS,
  __const int COLS
)
{
  int bx = get_group_id( 0), by = get_group_id( 1);
  int tx = get_local_id( 0), ty = get_local_id( 1);
  int A_Begin = by * BLOCK_SIZE * COLS;
  int B_Begin = bx * BLOCK_SIZE * COLS;
  int A_End = A_Begin + COLS - 1, A_idx, B_idx, k, o;
  PRECISION tmp, a_b_sum = 0.0, a2_sum = 0.0, b2_sum = 0.0, a = 0.0, b = 0.0;

  for( A_idx = A_Begin, B_idx = B_Begin; A_idx <= A_End; A_idx += BLOCK_SIZE, B_idx += BLOCK_SIZE)
  {
    As[ ty * BLOCK_SIZE + tx] = INPUT_A[ A_idx + ty * COLS + tx];
    Bs[ tx * BLOCK_SIZE + ty] = INPUT_A[ B_idx + ty * COLS + tx];
    barrier( CLK_LOCAL_MEM_FENCE);
#pragma unroll
    for( k = 0; k < BLOCK_SIZE; k++)
    {
      a = As[ ty * BLOCK_SIZE + k];
      b = Bs[ k * BLOCK_SIZE + tx];
      a_b_sum += a * b;
      a2_sum  += a * a;
      b2_sum  += b * b;
    }
    barrier( CLK_LOCAL_MEM_FENCE);
  }

  o = by * BLOCK_SIZE * ROWS + ty * ROWS + bx * BLOCK_SIZE + tx;
  OUTPUT[ o] = 2 * a_b_sum * native_recip( a2_sum + b2_sum);
}

//! @brief kernel for calculating pairwise Cosine similarity coefficients
//! @param OUTPUT buffer to place resulting distance matrix
//! @param INPUT_A, INPUT_A the two input matrices
//! @param As, Bs the local memory allocations
//! @param ROWS, COLS, ROWS, COLS the dimensionality
__kernel void Cosine
(
  __global PRECISION* OUTPUT,
  __const __global PRECISION* INPUT_A,
  __local PRECISION* As,
  __local PRECISION* Bs,
  __const int ROWS,
  __const int COLS
)
{
  int bx = get_group_id( 0), by = get_group_id( 1);
  int tx = get_local_id( 0), ty = get_local_id( 1);
  int A_Begin = by * BLOCK_SIZE * COLS;
  int B_Begin = bx * BLOCK_SIZE * COLS;
  int A_End = A_Begin + COLS - 1, A_idx, B_idx, k, o;
  PRECISION tmp, a_b_sum = 0.0, a2_sum = 0.0, b2_sum = 0.0, a = 0.0, b = 0.0;

  for( A_idx = A_Begin, B_idx = B_Begin; A_idx <= A_End; A_idx += BLOCK_SIZE, B_idx += BLOCK_SIZE)
  {
    As[ ty * BLOCK_SIZE + tx] = INPUT_A[ A_idx + ty * COLS + tx];
    Bs[ tx * BLOCK_SIZE + ty] = INPUT_A[ B_idx + ty * COLS + tx];
    barrier( CLK_LOCAL_MEM_FENCE);
#pragma unroll
    for( k = 0; k < BLOCK_SIZE; k++)
    {
      a = As[ ty * BLOCK_SIZE + k];
      b = Bs[ k * BLOCK_SIZE + tx];
      a_b_sum += a * b;
      a2_sum  += a * a;
      b2_sum  += b * b;
    }
    barrier( CLK_LOCAL_MEM_FENCE);
  }

  o = by * BLOCK_SIZE * ROWS + ty * ROWS + bx * BLOCK_SIZE + tx;
  OUTPUT[ o] = a_b_sum * native_recip( native_sqrt( a2_sum * b2_sum));
}

//! @brief kernel for calculating pairwise manhattan similarity coefficients
//! @param OUTPUT buffer to place resulting distance matrix
//! @param INPUT_A, INPUT_A the two input matrices
//! @param As, Bs the local memory allocations
//! @param ROWS, COLS, ROWS, COLS the dimensionality
__kernel void Manhattan
(
  __global PRECISION* OUTPUT,
  __const __global PRECISION* INPUT_A,
  __local PRECISION* As,
  __local PRECISION* Bs,
  __const int ROWS,
  __const int COLS
)
{
  int bx = get_group_id( 0), by = get_group_id( 1);
  int tx = get_local_id( 0), ty = get_local_id( 1);
  int A_Begin = by * BLOCK_SIZE * COLS;
  int B_Begin = bx * BLOCK_SIZE * COLS;
  int A_End = A_Begin + COLS - 1, A_idx, B_idx, k, o;
  PRECISION tmp, s = 0.0;

  for( A_idx = A_Begin, B_idx = B_Begin; A_idx <= A_End; A_idx += BLOCK_SIZE, B_idx += BLOCK_SIZE)
  {
    As[ ty * BLOCK_SIZE + tx] = INPUT_A[ A_idx + ty * COLS + tx];
    Bs[ tx * BLOCK_SIZE + ty] = INPUT_A[ B_idx + ty * COLS + tx];
    barrier( CLK_LOCAL_MEM_FENCE);
#pragma unroll
    for( k = 0; k < BLOCK_SIZE; k++)
    {
      tmp = fabs( As[ ty * BLOCK_SIZE + k] - Bs[ k * BLOCK_SIZE + tx]);
      s += tmp;
    }
    barrier( CLK_LOCAL_MEM_FENCE);
  }

  o = by * BLOCK_SIZE * ROWS + ty * ROWS + bx * BLOCK_SIZE + tx;
  OUTPUT[ o] = native_recip( 1 + s);
}

__kernel void RMSD
(
  __global PRECISION* OUTPUT,
  __const __global PRECISION* INPUT,
  __local PRECISION* Ys,
  __local PRECISION* Xs,
  __const int ROWS,
  __const int COLS,
  __const int NUM_ATOMS
)
{
  int bx = get_group_id( 0), by = get_group_id( 1);
  int tx = get_local_id( 0), ty = get_local_id( 1);
  int block_size = get_local_size( 0);
  int yBegin = by * block_size * COLS;
  int xBegin = bx * block_size * COLS;
  int yEnd = yBegin + COLS -1, y, x, k, o;
  PRECISION tmp, s = 0;

  for( y = yBegin, x = xBegin; y <= yEnd; y += block_size, x += block_size)
  {
    Ys[ ty * block_size + tx] = INPUT[ y + ty * COLS + tx];
    Xs[ tx * block_size + ty] = INPUT[ x + ty * COLS + tx];
    // transpose Xs to avoid bank conflicts
    barrier( CLK_LOCAL_MEM_FENCE);
    for( k = 0; k < block_size; k++)
    {
      tmp = Ys[ ty * block_size + k] - Xs[ k * block_size + tx];
      s += tmp * tmp;
    }
    barrier( CLK_LOCAL_MEM_FENCE);
  }

  o = by * block_size * ROWS + ty * ROWS + bx * block_size + tx;
  OUTPUT[ o] = native_sqrt( s / ( PRECISION)( NUM_ATOMS));
}

//! @brief kernel for calculating pairwise euclidean similarity coefficients
//! @param OUTPUT buffer to place resulting distance matrix
//! @param INPUT_A, INPUT_A the two input matrices
//! @param As, Bs the local memory allocations
//! @param ROWS, COLS, ROWS, COLS the dimensionality
__kernel void EuclideanSlow
(
  __global PRECISION* OUTPUT,
  __const  uint ROW,
  __const __global PRECISION* INPUT_A,
  __local PRECISION* As,
  __const int ROWS,
  __const int COLS
)
{
  PRECISION tmp, s = 0;
  int index = get_global_id( 0);

  for( int i = 0; i < COLS; ++i)
  {
    As[ i] = INPUT_A[ ROW * COLS + i];
  }
  barrier( CLK_LOCAL_MEM_FENCE);

  for( int col = 0; col < COLS; ++col)
  {
    tmp = As[ col] - INPUT_A[ index * COLS + col];
    s += tmp * tmp;
  }
  barrier( CLK_LOCAL_MEM_FENCE);

  OUTPUT[ index] = native_recip( 1 + native_sqrt( s));
}

//! @brief kernel for calculating pairwise tanimoto similarity coefficients
//! @param OUTPUT buffer to place resulting distance matrix
//! @param INPUT_A, INPUT_A the two input matrices
//! @param As, Bs the local memory allocations
//! @param ROWS, COLS, ROWS, COLS the dimensionality
__kernel void TanimotoSlow
(
  __global PRECISION* OUTPUT,
  __const  uint ROW,
  __const __global PRECISION* INPUT_A,
  __local PRECISION* As,
  __const int ROWS,
  __const int COLS
)
{
  PRECISION a_b_sum = 0, a2_sum = 0, b2_sum = 0, a, b;
  int index = get_global_id( 0);

  for( int i = 0; i < COLS; ++i)
  {
    As[ i] = INPUT_A[ ROW * COLS + i];
  }
  barrier( CLK_LOCAL_MEM_FENCE);

  for( int col = 0; col < COLS; ++col)
  {
    a = As[ col];
    b = INPUT_A[ index * COLS + col];
    a_b_sum += a * b;
    a2_sum  += a * a;
    b2_sum  += b * b;
  }
  barrier( CLK_LOCAL_MEM_FENCE);

  OUTPUT[ index] = a_b_sum * native_recip( a2_sum + b2_sum - a_b_sum);
}

//! @brief kernel for calculating pairwise Dice similarity coefficients
//! @param OUTPUT buffer to place resulting distance matrix
//! @param INPUT_A, INPUT_A the two input matrices
//! @param As, Bs the local memory allocations
//! @param ROWS, COLS, ROWS, COLS the dimensionality
__kernel void DiceSlow
(
  __global PRECISION* OUTPUT,
  __const  uint ROW,
  __const __global PRECISION* INPUT_A,
  __local PRECISION* As,
  __const int ROWS,
  __const int COLS
)
{
  PRECISION tmp, a_b_sum = 0.0, a2_sum = 0.0, b2_sum = 0.0, a = 0.0, b = 0.0;

  int index = get_global_id( 0);

  for( int i = 0; i < COLS; ++i)
  {
    As[ i] = INPUT_A[ ROW * COLS + i];
  }
  barrier( CLK_LOCAL_MEM_FENCE);

  for( int col = 0; col < COLS; ++col)
  {
    a = As[ col];
    b = INPUT_A[ index * COLS + col];
    a_b_sum += a * b;
    a2_sum  += a * a;
    b2_sum  += b * b;
  }
  barrier( CLK_LOCAL_MEM_FENCE);

  OUTPUT[ index] = 2 * a_b_sum * native_recip( a2_sum + b2_sum);
}

//! @brief kernel for calculating pairwise Cosine similarity coefficients
//! @param OUTPUT buffer to place resulting distance matrix
//! @param INPUT_A, INPUT_A the two input matrices
//! @param As, Bs the local memory allocations
//! @param ROWS, COLS, ROWS, COLS the dimensionality
__kernel void CosineSlow
(
  __global PRECISION* OUTPUT,
  __const  uint ROW,
  __const __global PRECISION* INPUT_A,
  __local PRECISION* As,
  __const int ROWS,
  __const int COLS
)
{
  PRECISION tmp, a_b_sum = 0.0, a2_sum = 0.0, b2_sum = 0.0, a = 0.0, b = 0.0;

  int index = get_global_id( 0);

  for( int i = 0; i < COLS; ++i)
  {
    As[ i] = INPUT_A[ ROW * COLS + i];
  }
  barrier( CLK_LOCAL_MEM_FENCE);

  for( int col = 0; col < COLS; ++col)
  {
    a = As[ col];
    b = INPUT_A[ index * COLS + col];
    a_b_sum += a * b;
    a2_sum  += a * a;
    b2_sum  += b * b;
  }
  barrier( CLK_LOCAL_MEM_FENCE);

  OUTPUT[ index] = a_b_sum * native_recip( native_sqrt( a2_sum * b2_sum));
}

//! @brief kernel for calculating pairwise manhattan similarity coefficients
//! @param OUTPUT buffer to place resulting distance matrix
//! @param INPUT_A, INPUT_A the two input matrices
//! @param As, Bs the local memory allocations
//! @param ROWS, COLS, ROWS, COLS the dimensionality
__kernel void ManhattanSlow
(
  __global PRECISION* OUTPUT,
  __const  uint ROW,
  __const __global PRECISION* INPUT_A,
  __local PRECISION* As,
  __const int ROWS,
  __const int COLS
)
{
  PRECISION tmp, s = 0.0;

  int index = get_global_id( 0);

  for( int i = 0; i < COLS; ++i)
  {
    As[ i] = INPUT_A[ ROW * COLS + i];
  }
  barrier( CLK_LOCAL_MEM_FENCE);

  for( int col = 0; col < COLS; ++col)
  {
    PRECISION b = INPUT_A[ index * COLS + col];
    tmp = fabs( As[ col] - b);
    s += tmp;
  }
  barrier( CLK_LOCAL_MEM_FENCE);

  OUTPUT[ index] = native_recip( 1 + s);
}

//! @brief kernel for calculating pairwise manhattan similarity coefficients
//! @param OUTPUT buffer to place resulting sums
//! @param INPUT_A, INPUT_A the two input matrices
//! @param As, Bs the local memory allocations
//! @param ROWS, COLS, ROWS, COLS the dimensionality
__kernel void CalcAllSums
(
  __global PRECISION* OUTPUT_A_B_SUM,
  __global PRECISION* OUTPUT_A2_SUM,
  __global PRECISION* OUTPUT_B2_SUM,
  __global PRECISION* OUTPUT_A_MINUS_B_SUM,
  __global PRECISION* OUTPUT_A_MINUS_B2_SUM,
  __const __global PRECISION* INPUT_A,
  __local PRECISION* As,
  __local PRECISION* Bs,
  __const int ROWS,
  __const int COLS
)
{
  int bx = get_group_id( 0), by = get_group_id( 1);
  int tx = get_local_id( 0), ty = get_local_id( 1);
  int A_Begin = by * BLOCK_SIZE * COLS;
  int B_Begin = bx * BLOCK_SIZE * COLS;
  int A_End = A_Begin + COLS - 1, A_idx, B_idx, k, o;
  PRECISION a, b, a_b_sum = 0.0, a2_sum = 0.0, b2_sum = 0.0, a_minus_b_sum = 0.0, a_minus_b2_sum = 0.0;

  for( A_idx = A_Begin, B_idx = B_Begin; A_idx <= A_End; A_idx += BLOCK_SIZE, B_idx += BLOCK_SIZE)
  {
    As[ ty * BLOCK_SIZE + tx] = INPUT_A[ A_idx + ty * COLS + tx];
    Bs[ tx * BLOCK_SIZE + ty] = INPUT_A[ B_idx + ty * COLS + tx];
    barrier( CLK_LOCAL_MEM_FENCE);
#pragma unroll
    for( k = 0; k < BLOCK_SIZE; k++)
    {
      a = As[ ty * BLOCK_SIZE + k];
      b = Bs[ k * BLOCK_SIZE + tx];

      a_b_sum        += a * b;
      a2_sum         += a * a;
      b2_sum         += b * b;
      a_minus_b_sum  += a - b;
      a_minus_b2_sum += a_minus_b_sum * a_minus_b_sum;
    }
    barrier( CLK_LOCAL_MEM_FENCE);
  }

  o = by * BLOCK_SIZE * ROWS + ty * ROWS + bx * BLOCK_SIZE + tx;

  OUTPUT_A_B_SUM[ o]        = a_b_sum;
  OUTPUT_A2_SUM[ o]         = a2_sum;
  OUTPUT_B2_SUM[ o]         = b2_sum;
  OUTPUT_A_MINUS_B_SUM[ o]  = a_minus_b_sum;
  OUTPUT_A_MINUS_B2_SUM[ o] = a_minus_b2_sum;
}

//! @brief kernel for calculating pairwise manhattan similarity coefficients
//! @param OUTPUT buffer to place resulting sums
//! @param INPUT_A, INPUT_A the two input matrices
//! @param As, Bs the local memory allocations
//! @param ROWS, COLS, ROWS, COLS the dimensionality
__kernel void DiceFromSums
(
          __global PRECISION* OUTPUT,
  __const __global PRECISION* INPUT_A_B_SUM,
  __const __global PRECISION* INPUT_A2_SUM,
  __const int ROWS,
  __const int COLS
)
{
  size_t index = get_global_id( 0);
  if( index >= ROWS) return;
  PRECISION a_b_sum, a2_sum, b2_sum;

  for( size_t i = 0; i < COLS; ++i)
  {
    a_b_sum = INPUT_A_B_SUM[ index * COLS + i];
    a2_sum  = INPUT_A2_SUM[ index * COLS + i];
    b2_sum  = INPUT_A2_SUM[ index * COLS + i];

    OUTPUT[ index * COLS + i] = 2 * a_b_sum * native_recip( a2_sum + b2_sum);
  }
}

//! @brief kernel for calculating pairwise manhattan similarity coefficients
//! @param OUTPUT buffer to place resulting sums
//! @param INPUT_A, INPUT_A the two input matrices
//! @param As, Bs the local memory allocations
//! @param ROWS, COLS, ROWS, COLS the dimensionality
__kernel void CosineFromSums
(
          __global PRECISION* OUTPUT,
  __const __global PRECISION* INPUT_A_B_SUM,
  __const __global PRECISION* INPUT_A2_SUM,
  __const int ROWS,
  __const int COLS
)
{
  size_t index = get_global_id( 0);
  if( index >= ROWS) return;
  PRECISION a_b_sum, a2_sum, b2_sum;

  for( size_t i = 0; i < COLS; ++i)
  {
    a_b_sum = INPUT_A_B_SUM[ index * COLS + i];
    a2_sum  = INPUT_A2_SUM[ index * COLS + i];
    b2_sum  = INPUT_A2_SUM[ index * COLS + i];

    OUTPUT[ index * COLS + i] = a_b_sum * native_recip( native_sqrt( a2_sum * b2_sum));
  }
}

//! @brief kernel for calculating pairwise manhattan similarity coefficients
//! @param OUTPUT buffer to place resulting sums
//! @param INPUT_A, INPUT_A the two input matrices
//! @param As, Bs the local memory allocations
//! @param ROWS, COLS, ROWS, COLS the dimensionality
__kernel void EuclideanFromSums
(
          __global PRECISION* OUTPUT,
  __const __global PRECISION* INPUT_A_MINUS_B2_SUM,
  __const int ROWS,
  __const int COLS
)
{
  size_t index = get_global_id( 0);
  if( index >= ROWS) return;
  PRECISION a_minus_b2_sum;

  for( size_t i = 0; i < COLS; ++i)
  {
    a_minus_b2_sum = INPUT_A_MINUS_B2_SUM[ index * COLS + i];

    OUTPUT[ index * COLS + i] = native_recip( 1 + native_sqrt( a_minus_b2_sum));
  }
}

//! @brief kernel for calculating pairwise manhattan similarity coefficients
//! @param OUTPUT buffer to place resulting sums
//! @param INPUT_A, INPUT_A the two input matrices
//! @param As, Bs the local memory allocations
//! @param ROWS, COLS, ROWS, COLS the dimensionality
__kernel void ManhattanFromSums
(
          __global PRECISION* OUTPUT,
  __const __global PRECISION* INPUT_A_MINUS_B_SUM,
  __const int ROWS,
  __const int COLS
)
{
  size_t index = get_global_id( 0);
  if( index >= ROWS) return;
  PRECISION a_minus_b_sum;

  for( size_t i = 0; i < COLS; ++i)
  {
    a_minus_b_sum = INPUT_A_MINUS_B_SUM[ index * COLS + i];

    OUTPUT[ index * COLS + i] = native_recip( 1 + fabs( a_minus_b_sum));
  }
}

//! @brief kernel for calculating pairwise manhattan similarity coefficients
//! @param OUTPUT buffer to place resulting sums
//! @param INPUT_A, INPUT_A the two input matrices
//! @param As, Bs the local memory allocations
//! @param ROWS, COLS, ROWS, COLS the dimensionality
__kernel void TanimotoFromSums
(
          __global PRECISION* OUTPUT,
  __const __global PRECISION* INPUT_A_B_SUM,
  __const __global PRECISION* INPUT_A2_SUM,
  __const int ROWS,
  __const int COLS
)
{
  size_t index = get_global_id( 0);
  if( index >= ROWS) return;
  PRECISION a_b_sum, a2_sum, b2_sum;

  for( size_t i = 0; i < COLS; ++i)
  {
    a_b_sum = INPUT_A_B_SUM[ index * COLS + i];
    a2_sum  = INPUT_A2_SUM[ index * COLS + i];
    b2_sum  = INPUT_A2_SUM[ index * COLS + i];

    OUTPUT[ index * COLS + i] = a_b_sum * native_recip( a2_sum + b2_sum - a_b_sum);
  }
}
