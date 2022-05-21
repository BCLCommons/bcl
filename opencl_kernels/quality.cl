// IMPORTANT: AMD currently has a problem compiling kernels with trigonometric functions in them
// If the BCL fails to start on your machine and you have either an AMD card installed, or have used
// -opencl AMD_Accelerated_Parallel_Processing as a flag, then please delete this file.  The BCL should run properly
// using CPU versions of the code in question.  Alternatively, install the NVIDIA drivers.  Alternatively, just remove
// the RMSDFromCovarianceMatrix function below

// this is only included, so that eclipse c++ editor can be used and parses information correctly
// this will not be interpreted, when the program is compiled with an OpenCL compiler
#ifndef __OPENCL_VERSION__
 #include "kernel.h"
#endif

//! @brief calculate the the rmsd between two sets of coordinates
//! @param COORDIANTES the coordinates that are usually superimposed as matrix with local_size(0) cols and NUMBER_COORDINATES rows
//! @param REFERENCE_COORDINATES the coordinates that are the reference as matrix with local_size(0) cols and NUMBER_COORDINATES rows
//! @param NUMBER_COORDINATES the number of coordinates (rows) in matrices
//! @param RMSD output root mean square deviation of row coordinates
//! @param SHARED accumulator for reduction
__kernel void RealSpaceRMSD
(
  const __global PRECISION* COORDINATES,
  const __global PRECISION* REFERENCE_COORDINATES,
  const          uint       NUMBER_COORDINATES,
        __global PRECISION* RMSD,
        __local  PRECISION* SHARED
)
{
  // local id
  const size_t local_index = get_local_id( 0);

  // block size
  const size_t block_size = get_local_size( 0);

  // accumulator
  PRECISION accumulator = 0.0;

  for( size_t index = local_index; index < NUMBER_COORDINATES; index += block_size)
  {
    const PRECISION diff = COORDINATES[ index] - REFERENCE_COORDINATES[ index];
    accumulator += diff * diff;
  }

  SHARED[ local_index] = accumulator;

  // Synchronize to make sure that shared mem if full
  barrier( CLK_LOCAL_MEM_FENCE);

  // parallel reduction
  for( size_t offset = block_size / 2; offset > 0; offset = offset / 2)
  {
    if( local_index < offset)
    {
      SHARED[ local_index] += SHARED[ local_index + offset];
    }
    barrier( CLK_LOCAL_MEM_FENCE);
  }

  if( local_index == 0)
  {
    RMSD[ 0] = SHARED[ 0];
  }
}

//! @brief calculate the average position of a set of coordinates
//! @param COORDINATES coordinates to consider
//! @param NUM_COORDINATES number of rows in array of coordinates
//! @param CENTER the calculated center
__kernel void CoordinatesCenter
(
  const __global PRECISION* COORDINATES,
  const          uint       NUM_COORDINATES,
        __global PRECISION* CENTER,
        __local  PRECISION* SHARED
)
{
  // Perform parallel reduction
  __const size_t local_index_x = get_local_id( 0);
  __const size_t local_index_y = get_local_id( 1);

  __const size_t block_size = get_global_size( 1);

  __const size_t step = block_size * block_size;

  __const size_t max = block_size * NUM_COORDINATES;

  __const size_t start = local_index_y * block_size + local_index_x;
  PRECISION accumulator = 0.0;

  // Loop sequentially over chunks of input vector
  for( size_t coord_index = start; coord_index < max; coord_index += step)
  {
    accumulator += COORDINATES[ coord_index];
  }

  SHARED[ start] = accumulator;
  barrier( CLK_LOCAL_MEM_FENCE);

  for( size_t offset = block_size / 2; offset > 0; offset = offset / 2)
  {
    if( local_index_y < offset)
    {
      SHARED[ local_index_y * block_size + local_index_x] += SHARED[ ( local_index_y + offset) * block_size + local_index_x];
    }
    barrier( CLK_LOCAL_MEM_FENCE);
  }
  if( local_index_y == 0)
  {
    CENTER[ local_index_x] = SHARED[ local_index_x] / NUM_COORDINATES;
  }
}

//! @brief calculate the centers of all fragments starting from row (offset)
//! global id 0 is the coordinate component (x, y, z)
//! global id 1 is the fragment number that is considered
//! @param COORDINATES the coordinate matrix
//! @param OFFSET first row of first fragment
//! @param FRAGMENT_LENGTH number of coordinate rows per fragment
//! @param CENTERS output
__kernel void FragmentCenters
(
  const __global PRECISION* COORDINATES,
  const          uint       FRAGMENT_LENGTH,
  const          uint       NUMBER_FRAGMENTS,
        __global PRECISION* CENTERS
)
{
  // index of coordinate x, y, z
  const size_t thread_x        = get_global_id( 0);
  const size_t number_cols     = get_global_size( 0);
  const size_t fragment_number = get_global_id( 1);

  if( fragment_number >= NUMBER_FRAGMENTS)
  {
    return;
  }

  // row to calculate center over
  const size_t start_row = fragment_number;
  const size_t end_row   = fragment_number + FRAGMENT_LENGTH;

  // accumulator
  PRECISION coord_sum = 0.0;

  // iterate over rows for current fragment (given by fragment_number
  for( size_t row = start_row; row < end_row; ++ row)
  {
    coord_sum += COORDINATES[ row * number_cols + thread_x];
  }

  CENTERS[ start_row * number_cols + thread_x] = coord_sum / FRAGMENT_LENGTH;
}

//! @brief calculate the average position of a set of coordinates according to a selection
//! @param COORDINATES coordinates to consider
//! @param SELECTION for each vector in COORDINATES if it should be considered
//! @param NUM_COORDINATES number of rows in array of coordinates
//! @param NUM_SELECTION number of selected coordinates
//! @param CENTER the calculated center
__kernel void CoordinatesCenterSelections
(
  const __global PRECISION* COORDINATES,
  const __global int*       SELECTIONS,
  const          uint       NUM_COLS,
  const          uint       NUM_COORDINATES,
        __global PRECISION* CENTERS
)
{
  // Perform parallel reduction
  const size_t index_x         = get_local_id( 0);
  const size_t fragment_number = get_global_id( 1);
  const size_t num_cols_select = get_global_size( 1);

  PRECISION accumulator = 0.0;
  size_t counter = 0;

  // Loop sequentially over input coordinates
  for( size_t row = 0; row < NUM_COORDINATES; ++row)
  {
    // consider coordinate if selected
    if( SELECTIONS[ row * num_cols_select + fragment_number] > 0)
    {
      accumulator += COORDINATES[ row * NUM_COLS + index_x];
      ++counter;
    }
  }

  CENTERS[ fragment_number * NUM_COLS + index_x] = accumulator / counter;
}

//! @brief compute covariance matrix of two sets of coordinates COORDINATES on REFERENCE_COORDINATES
//! both coordinate sets are translated to the center of mass
//! @param COORDINATES set of coordinates
//! @param REFERENCE_COORDINATES set of coordinates
//! @param NUM_COORDINATES number of rows in array of coordinates
//! @param CENTER_A the center of COORDINATES
//! @param CENTER_B the center of REFERENCE_COORDINATES
//! @param COVERIANCE_MATRIX resulting covariance matrix
//! @param SQUARE_NORM_CENTERED_COORDINATES optional pointer to which the square norm of the centered coordinates a will be depsosited
//! @param SQUARE_NORM_CENTERED_REFERENCE_COORDINATES optional pointer to which the square norm of the centered coordinates b will be depsosited
__kernel void BuildCovarianceMatrix
(
  const __global PRECISION* COORDINATES,
  const __global PRECISION* REFERENCE_COORDINATES,
  const          uint       NUM_COORDINATES,
  const          uint       NUM_COLS,
  const __global PRECISION* CENTER_A,
  const __global PRECISION* CENTER_B,
        __global PRECISION* COVERIANCE_MATRIX,
        __global PRECISION* SQUARE_NORM_CENTERED_COORDINATES,
        __global PRECISION* SQUARE_NORM_CENTERED_REFERENCE_COORDINATES,
        __local  PRECISION* As,
        __local  PRECISION* Bs,
        __local  PRECISION* SH_CENTER_A,
        __local  PRECISION* SH_CENTER_B,
        __local  PRECISION* INNER_PRODUCT_A,
        __local  PRECISION* INNER_PRODUCT_B
)
{
  // Thread index
  __const size_t tx = get_local_id( 0);
  __const size_t ty = get_local_id( 1);

  __const size_t block_size = get_local_size( 0);

  INNER_PRODUCT_A[ ty * block_size + tx] = 0.0;
  INNER_PRODUCT_B[ ty * block_size + tx] = 0.0;

  // copy center to shared memory
  if( ty == 0)
  {
    SH_CENTER_A[ tx] = CENTER_A[ tx];
    SH_CENTER_B[ tx] = CENTER_B[ tx];
  }
  // Synchronize to make sure the centers are loaded
  barrier( CLK_LOCAL_MEM_FENCE);

  // Csub is used to store the element of the block sub-matrix
  // that is computed by the thread
  PRECISION c_sub = 0.0;

  const size_t num_rows = ( ( ( NUM_COORDINATES - 1) / block_size) + 1) * block_size;

  // Loop over all the sub-matrices of A and B
  // required to compute the block sub-matrix
  for( size_t row = 0; row < num_rows; row += block_size)
  {
    // Load the matrices from device memory
    // to shared memory; each thread loads
    // one element of each matrix
    if( row + ty < NUM_COORDINATES)
    {
      As[ ty * NUM_COLS + tx] = COORDINATES[ ( row + ty) * block_size + tx] - SH_CENTER_A[ tx];
      Bs[ ty * NUM_COLS + tx] = REFERENCE_COORDINATES[ ( row + ty) * block_size + tx] - SH_CENTER_B[ tx];
    }
    else
    {
      As[ ty * NUM_COLS + tx] = 0.0;
      Bs[ ty * NUM_COLS + tx] = 0.0;
    }
    INNER_PRODUCT_A[ ty * block_size + tx] += pown( As[ ty * NUM_COLS + tx], 2);
    INNER_PRODUCT_B[ ty * block_size + tx] += pown( Bs[ ty * NUM_COLS + tx], 2);

    // Synchronize to make sure the matrices are loaded
    barrier( CLK_LOCAL_MEM_FENCE);

    // Multiply the two matrices together;
    // each thread computes one element
    // of the block sub-matrix
//    #pragma unroll
    for( size_t k = 0; k < block_size; ++k)
    {
      c_sub += As[ k * NUM_COLS + ty] * Bs[ k * NUM_COLS + tx];
    }

    // Synchronize to make sure that the preceding
    // computation is done before loading two new
    // sub-matrices of A and B in the next iteration
    barrier( CLK_LOCAL_MEM_FENCE);
  }

  // Write the block sub-matrix to device memory;
  // each thread writes one element
  COVERIANCE_MATRIX[ ty * block_size + tx] = c_sub;

  const size_t local_index = ty * block_size + tx;

  // parallel reduction of the inner product
  for( size_t offset = block_size * block_size / 2; offset > 0; offset = offset / 2)
  {
    if( local_index < offset)
    {
      INNER_PRODUCT_A[ local_index] += INNER_PRODUCT_A[ local_index + offset];
      INNER_PRODUCT_B[ local_index] += INNER_PRODUCT_B[ local_index + offset];
    }
    barrier( CLK_LOCAL_MEM_FENCE);
  }
  if( local_index == 0)
  {
    SQUARE_NORM_CENTERED_COORDINATES[ 0] = INNER_PRODUCT_A[ 0];
    SQUARE_NORM_CENTERED_REFERENCE_COORDINATES[ 0] = INNER_PRODUCT_B[ 0];
  }
}

//! @brief compute covariance matrix of two sets of coordinates COORDINATES on REFERENCE_COORDINATES
//! both coordinate sets are translated to the center of mass
//! @param COORDINATES set of coordinates
//! @param REFERENCE_COORDINATES set of coordinates
//! @param NUM_COORDINATES number of rows in array of coordinates
//! @param CENTER_COORDS the center of COORDINATES
//! @param CENTER_REF_COORDS the center of REFERENCE_COORDINATES
//! @param COVERIANCE_MATRIX resulting covariance matrix
//! @param SQUARE_NORM_CENTERED_COORDINATES optional pointer to which the square norm of the centered coordinates a will be depsosited
//! @param SQUARE_NORM_CENTERED_REFERENCE_COORDINATES optional pointer to which the square norm of the centered coordinates b will be depsosited
__kernel void BuildCovarianceMatrixFragments
(
  const __global PRECISION* COORDINATES,
  const __global PRECISION* REFERENCE_COORDINATES,
  const __global PRECISION* CENTER_COORDS,
  const __global PRECISION* CENTER_REF_COORDS,
  const          uint       NUM_COLS,
  const          uint       FRAGMENT_LENGTH,
        __global PRECISION* COVERIANCE_MATRIX,
  const          uint       MATRIX_NUM_ROWS,
        __global PRECISION* SQUARE_NORM_CENTERED_COORDINATES,
        __global PRECISION* SQUARE_NORM_CENTERED_REFERENCE_COORDINATES,
        __local  PRECISION* As,
        __local  PRECISION* Bs,
        __local  PRECISION* SH_CENTER_A,
        __local  PRECISION* SH_CENTER_B,
        __local  PRECISION* INNER_PRODUCT_A,
        __local  PRECISION* INNER_PRODUCT_B
)
{
  // Thread index
  const size_t tx = get_local_id( 0);
  const size_t ty = get_local_id( 1);
  const size_t fragment_number = get_group_id( 2);

  const size_t block_size_x = get_local_size( 0);
  const size_t block_size_y = get_local_size( 1);

  INNER_PRODUCT_A[ ty * block_size_x + tx] = 0.0;
  INNER_PRODUCT_B[ ty * block_size_x + tx] = 0.0;

  // copy center to shared memory
  if( ty == 0)
  {
    SH_CENTER_A[ tx] = CENTER_COORDS[ fragment_number * NUM_COLS + tx];
    SH_CENTER_B[ tx] = CENTER_REF_COORDS[ fragment_number * NUM_COLS + tx];
  }
  // Synchronize to make sure the centers are loaded
  barrier( CLK_LOCAL_MEM_FENCE);

  // Csub is used to store the element of the block sub-matrix
  // that is computed by the thread
  PRECISION c_sub = 0.0;

  const size_t start_row = fragment_number;
  const size_t num_rows = ( ( ( FRAGMENT_LENGTH - 1) / block_size_y) + 1) * block_size_y;

  // Loop over all the sub-matrices of A and B
  // required to compute the block sub-matrix
  for( size_t row = 0; row < num_rows; row += block_size_y)
  {
    // Load the matrices from device memory
    // to shared memory; each thread loads
    // one element of each matrix
    if( row + ty < FRAGMENT_LENGTH)
    {
      As[ ty * NUM_COLS + tx] = COORDINATES[ ( start_row + row + ty) * NUM_COLS + tx] - SH_CENTER_A[ tx];
      Bs[ ty * NUM_COLS + tx] = REFERENCE_COORDINATES[ ( start_row + row + ty) * NUM_COLS + tx] - SH_CENTER_B[ tx];
    }
    else
    {
      As[ ty * NUM_COLS + tx] = 0.0;
      Bs[ ty * NUM_COLS + tx] = 0.0;
    }
    INNER_PRODUCT_A[ ty * block_size_x + tx] += pown( As[ ty * NUM_COLS + tx], 2);
    INNER_PRODUCT_B[ ty * block_size_x + tx] += pown( Bs[ ty * NUM_COLS + tx], 2);

    // Synchronize to make sure the matrices are loaded
    barrier( CLK_LOCAL_MEM_FENCE);

    // Multiply the two matrices together;
    // each thread computes one element
    // of the block sub-matrix
//    #pragma unroll
    for( size_t k = 0; k < block_size_x; ++k)
    {
      c_sub += As[ k * NUM_COLS + ty] * Bs[ k * NUM_COLS + tx];
    }

    // Synchronize to make sure that the preceding
    // computation is done before loading two new
    // sub-matrices of A and B in the next iteration
    barrier( CLK_LOCAL_MEM_FENCE);
  }

  // Write the block sub-matrix to device memory;
  // each thread writes one element
  if( ty < MATRIX_NUM_ROWS)
  {
    COVERIANCE_MATRIX[ fragment_number * MATRIX_NUM_ROWS * NUM_COLS + ty * NUM_COLS + tx] = c_sub;
  }

  __const size_t local_index = ty * block_size_x + tx;

  // parallel reduction of the inner product
  for( size_t offset = block_size_x * block_size_x / 2; offset > 0; offset = offset / 2)
  {
    if( local_index < offset)
    {
      INNER_PRODUCT_A[ local_index] += INNER_PRODUCT_A[ local_index + offset];
      INNER_PRODUCT_B[ local_index] += INNER_PRODUCT_B[ local_index + offset];
    }
    barrier( CLK_LOCAL_MEM_FENCE);
  }
  if( local_index == 0)
  {
    SQUARE_NORM_CENTERED_COORDINATES[ fragment_number] = INNER_PRODUCT_A[ 0];
    SQUARE_NORM_CENTERED_REFERENCE_COORDINATES[ fragment_number] = INNER_PRODUCT_B[ 0];
  }
}

//! @brief compute covariance matrix of two sets of coordinates COORDINATES on REFERENCE_COORDINATES for a selection
//! both coordinate sets are translated to the center of mass
//! @param COORDINATES set of coordinates
//! @param REFERENCE_COORDINATES set of coordinates
//! @param SELECTION for each vector in COORDINATES if it should be considered
//! @param NUM_COORDINATES number of rows in array of coordinates
//! @param NUM_SELECTION number of selected coordinates
//! @param CENTER_A the center of COORDINATES
//! @param CENTER_B the center of REFERENCE_COORDINATES
//! @param COVERIANCE_MATRIX resulting covariance matrix
//! @param NUMBER_ROWS_MATRIX number of rows per covariance matrix
__kernel void BuildCovarianceMatrixSelections
(
  const __global PRECISION* COORDINATES,
  const __global PRECISION* REFERENCE_COORDINATES,
  const __global int*       SELECTIONS,
  const          uint       NUM_COLS,
  const          uint       NUM_COORDINATES,
  const          uint       NUM_COLS_SELECTION,
  const __global PRECISION* CENTER_A,
  const __global PRECISION* CENTER_B,
        __global PRECISION* COVERIANCE_MATRIX,
  const          uint       NUMBER_ROWS_MATRIX,
        __local  PRECISION* As,
        __local  PRECISION* Bs,
        __local  PRECISION* SH_CENTER_A,
        __local  PRECISION* SH_CENTER_B
)
{
  // Thread index
  const size_t tx              = get_local_id( 0);
  const size_t ty              = get_local_id( 1);
  const size_t block_size      = get_local_size( 0);
  const size_t fragment_number = get_group_id( 2);

  // copy center to shared memory
  if( ty == 0)
  {
    SH_CENTER_A[ tx] = CENTER_A[ fragment_number * NUM_COLS + tx];
    SH_CENTER_B[ tx] = CENTER_B[ fragment_number * NUM_COLS + tx];
  }
  // Synchronize to make sure the centers are loaded
  barrier( CLK_LOCAL_MEM_FENCE);

  // Csub is used to store the element of the block sub-matrix
  // that is computed by the thread
  PRECISION c_sub = 0.0;

  const size_t num_rows = ( ( ( NUM_COORDINATES - 1) / block_size) + 1) * block_size;

  // Loop over all the sub-matrices of A and B
  // required to compute the block sub-matrix
  for( size_t row = 0; row < num_rows; row += block_size)
  {
    // Load the matrices from device memory
    // to shared memory; each thread loads
    // one element of each matrix
    // ignore everything beyond the number of coordinates and items that are ignored
    if( row + ty < NUM_COORDINATES && SELECTIONS[ ( row + ty) * NUM_COLS_SELECTION + fragment_number] > 0)
    {
      As[ ty * NUM_COLS + tx] = COORDINATES[ ( row + ty) * NUM_COLS + tx] - SH_CENTER_A[ tx];
      Bs[ ty * NUM_COLS + tx] = REFERENCE_COORDINATES[ ( row + ty) * NUM_COLS + tx] - SH_CENTER_B[ tx];
    }
    else
    {
      As[ ty * NUM_COLS + tx] = 0.0;
      Bs[ ty * NUM_COLS + tx] = 0.0;
    }

    // Synchronize to make sure the matrices are loaded
    barrier( CLK_LOCAL_MEM_FENCE);

    // Multiply the two matrices together;
    // each thread computes one element
    // of the block sub-matrix
//    #pragma unroll
    for( size_t k = 0; k < block_size; ++k)
    {
      c_sub += As[ k * NUM_COLS + tx] * Bs[ k * NUM_COLS + ty];
    }

    // Synchronize to make sure that the preceding
    // computation is done before loading two new
    // sub-matrices of A and B in the next iteration
    barrier( CLK_LOCAL_MEM_FENCE);
  }

  // Write the block sub-matrix to device memory;
  // each thread writes one element
  if( ty < NUMBER_ROWS_MATRIX)
  {
    COVERIANCE_MATRIX[ fragment_number * NUMBER_ROWS_MATRIX * NUM_COLS + ty * NUM_COLS + tx] = c_sub;
  }
}

//! @brief identify all coordinates that are below a given cutoff given a transformation
//! @param COORDINATES set of coordinates
//! @param REFERENCE_COORDINATES set of coordinates
//! @param SELECTIONS for each vector in COORDINATES and fragment number if it should be considered
//! @param NUM_COLS number of cols in coordinates and transformation matrices
//! @param NUM_COORDINATES number of rows in array of coordinates
//! @param NUM_COLS_SELECTION number of cols in selection matrix
//! @param TRANSFORMATIONS transformation matrix for each selection
//! @param NUMBER_ROWS_MATRIX number of rows in transformation matrices
//! @param SQUARE_DISTANCE_CUTOFF the cutoff for selection
__kernel void CoordinateSelectionsBelowCutoff
(
  const __global PRECISION* COORDINATES,
  const __global PRECISION* REFERENCE_COORDINATES,
        __global int*       SELECTIONS,
  const          uint       NUM_COLS,
  const          uint       NUM_COORDINATES,
  const          uint       NUM_COLS_SELECTION,
  const __global PRECISION* TRANSFORMATIONS,
  const          uint       NUMBER_ROWS_MATRIX,
  const          PRECISION  SQUARE_DISTANCE_CUTOFF,
        __local  PRECISION* TRANS_SHARED
)
{
  const size_t copy_index      = get_local_id( 0);
  const size_t row             = get_global_id( 0);
  const size_t fragment_number = get_global_id( 1);

  // copy transformation
  TRANS_SHARED[ copy_index] = TRANSFORMATIONS[ fragment_number * NUMBER_ROWS_MATRIX * NUM_COLS + copy_index];
  barrier( CLK_LOCAL_MEM_FENCE);

  //
  if( row >= NUM_COORDINATES)
  {
    return;
  }

  PRECISION old_coord[ 3];

  // transform coordinates
  for( size_t i = 0; i < 3; ++i)
  {
    old_coord[ i] = COORDINATES[ row * NUM_COLS + i];
  }

  PRECISION dist = 0.0;
  for( size_t col = 0; col < 3; ++col)
  {
    PRECISION acc = TRANS_SHARED[ 3 * NUM_COLS + col];

    for( size_t j = 0; j < 3; ++j)
    {
      acc += old_coord[ j] * TRANS_SHARED[ j * NUM_COLS + col];
    }
    acc -= REFERENCE_COORDINATES[ row * NUM_COLS + col];
    dist += pown( acc, 2);
  }

  SELECTIONS[ row * NUM_COLS_SELECTION + fragment_number] = islessequal( dist, SQUARE_DISTANCE_CUTOFF);
}

//! @brief using a set of covariance matrices, calculate the RMSD for each matrix/fragment
//! @param COVERIANCE_MATRIX covarinace matrices for each fragment one
//! @param FRAGMENT_LENGTH number of coordinates per fragment
//! @param NUMBER_FRAGMENTS number of fragments which is also the number of covarinace matrices
//! @param MATRIX_NUM_ROWS number of rows per covariance matrix
//! @param MATRIX_NUM_COLS number of cols per covariance matrix
//! @param SQUARE_NORM_CENTERED_COORDINATES sum of square norms
//! @param SQUARE_NORM_CENTERED_REFERENCE_COORDINATES sum of square norms
//! @param RMSDS output of rmsds vector - one element per fragment
__kernel void RMSDFromCovarianceMatrix
(
  const __global PRECISION* COVERIANCE_MATRIX,
  const               uint  FRAGMENT_LENGTH,
  const               uint  NUMBER_FRAGMENTS,
  const               uint  MATRIX_NUM_ROWS,
  const               uint  MATRIX_NUM_COLS,
  const __global PRECISION* SQUARE_NORM_CENTERED_COORDINATES,
  const __global PRECISION* SQUARE_NORM_CENTERED_REFERENCE_COORDINATES,
        __global PRECISION* RMSDS
)
{
  const size_t fragment_number = get_global_id( 0);

  if( fragment_number >= NUMBER_FRAGMENTS)
  {
    return;
  }

  const __global PRECISION *current_matrix = COVERIANCE_MATRIX + fragment_number * MATRIX_NUM_ROWS * MATRIX_NUM_COLS;

  const PRECISION o_00a = current_matrix[ 0 * MATRIX_NUM_COLS + 0];
  const PRECISION o_01b = current_matrix[ 0 * MATRIX_NUM_COLS + 1];
  const PRECISION o_02c = current_matrix[ 0 * MATRIX_NUM_COLS + 2];
  const PRECISION o_10d = current_matrix[ 1 * MATRIX_NUM_COLS + 0];
  const PRECISION o_11e = current_matrix[ 1 * MATRIX_NUM_COLS + 1];
  const PRECISION o_12f = current_matrix[ 1 * MATRIX_NUM_COLS + 2];
  const PRECISION o_20g = current_matrix[ 2 * MATRIX_NUM_COLS + 0];
  const PRECISION o_21h = current_matrix[ 2 * MATRIX_NUM_COLS + 1];
  const PRECISION o_22i = current_matrix[ 2 * MATRIX_NUM_COLS + 2];

  // calculate determinante
  const PRECISION determinant = o_00a * o_11e * o_22i
                              + o_01b * o_12f * o_20g
                              + o_02c * o_10d * o_21h
                              - o_02c * o_11e * o_20g
                              - o_01b * o_10d * o_22i
                              - o_00a * o_12f * o_21h;

  const int chi = determinant < 1E-10 ? -1 : 1;

  // multiply with its transposed to symmetrize
  const PRECISION m_00a = o_00a * o_00a + o_01b * o_01b + o_02c * o_02c;
  const PRECISION m_01b = o_00a * o_10d + o_01b * o_11e + o_02c * o_12f;
  const PRECISION m_02c = o_00a * o_20g + o_01b * o_21h + o_02c * o_22i;
  const PRECISION m_10d = m_01b;
  const PRECISION m_11e = o_10d * o_10d + o_11e * o_11e + o_12f * o_12f;
  const PRECISION m_12f = o_10d * o_20g + o_11e * o_21h + o_12f * o_22i;
  const PRECISION m_20g = m_02c;
  const PRECISION m_21h = m_12f;
  const PRECISION m_22i = o_20g * o_20g + o_21h * o_21h + o_22i * o_22i;

  // calculate eigenvalues
  // matrix a b c \n d e f \n g h i
  // characteristic polynomial
  // a = a + e + i
  // b =
  //Use the equation from above to get your cubic equation and combine all constant terms possible to
  //give you a reduced equation we will use a, b, c and d to denote the coefficients of this equation.
  //Eqn = a*lambda^3 + b*lambda^2 + c*lambda + d = 0
  const PRECISION a = -1.0;
  const PRECISION b = m_00a + m_11e + m_22i; // diagonal
  const PRECISION c = m_01b * m_10d + m_02c * m_20g +
                      m_12f * m_21h - m_00a * m_11e -
                      m_00a * m_22i - m_11e * m_22i;
  const PRECISION d = m_00a * m_11e * m_22i
                    + m_01b * m_12f * m_20g
                    + m_02c * m_10d * m_21h
                    - m_02c * m_11e * m_20g
                    - m_01b * m_10d * m_22i
                    - m_00a * m_12f * m_21h;

  const PRECISION x = c / a - pown( b, 2) / ( 3 * pown( a, 2));
  const PRECISION y = ( 2 * pown( b, 3) / pown( a, 3) - 9 * b * c / pown( a, 2) + 27 * d / a) / 27;
  const PRECISION z = pown( y, 2) / 4 + pown( x, 3) / 27;

  const PRECISION i = sqrt( pown( y, 2) / 4 - z);
  const PRECISION j = cbrt( i);
  const PRECISION k = acos( -y / ( 2 * i));
  const PRECISION m = cos( k / 3);
  const PRECISION n = sqrt( 3.0) * sin( k / 3);
  const PRECISION p = -b / ( 3 * a);

  // eigenvalues and sorting
  PRECISION eigenvalues[ 3] = { sqrt( 2 * j * m + p), sqrt( -j * ( m + n) + p), sqrt( -j * ( m - n) + p)};
  PRECISION tmp;
  if( eigenvalues[ 2] < eigenvalues[ 1])
  {
    tmp = eigenvalues[ 2];
    eigenvalues[ 2] = eigenvalues[ 1];
    eigenvalues[ 1] = tmp;
  }
  if( eigenvalues[ 1] < eigenvalues[ 0])
  {
    tmp = eigenvalues[ 1];
    eigenvalues[ 1] = eigenvalues[ 0];
    eigenvalues[ 0] = tmp;
  }
  if( eigenvalues[ 2] < eigenvalues[ 1])
  {
    tmp = eigenvalues[ 2];
    eigenvalues[ 2] = eigenvalues[ 1];
    eigenvalues[ 1] = tmp;
  }

  eigenvalues[ 0] *= chi;
  PRECISION square_deviation = 0.0;
  square_deviation += SQUARE_NORM_CENTERED_COORDINATES[ fragment_number];
  square_deviation += SQUARE_NORM_CENTERED_REFERENCE_COORDINATES[ fragment_number];
  square_deviation -= 2 * ( eigenvalues[ 0] + eigenvalues[ 1] + eigenvalues[ 2]);

  RMSDS[ fragment_number] = sqrt( max( square_deviation, 0.0) / FRAGMENT_LENGTH);
}

//! @brief update the best selections and transformations for each fragment, if the current selection has a higher count than the best selection
//! @param TRANSFORMATIONS current transformations for given current selections
//! @param SELECTIONS current selection collected with the transformations
//! @param TRANSFORMATIONS_BEST best transformation matrices so far
//! @param SELECTIONS_BEST best selections
//! @param NUM_COLS number of columns in transformation matrix
//! @param NUM_FRAGMENTS number of fragments equal to the number of transformation matrices
//! @param NUM_COORDINATES number of coordinates which the number of rows in the selections matrix
//! @param COUNTS_BEST sum of selections in best selections
//! @param DIFFERENCES difference of the current to the best selections sum - if negative, no update and no improvement, if positive update and improvement
//!                    this is used as termination criteria (i.e. no improvement) for GDT iterative superimposition algorithm
//! @param SHARED_SUM_SELECTION buffer for reduction
//! @param SHARED_SUM_SELECTION_BEST buffer for reduction
__kernel void UpdateSelections
(
  const __global PRECISION *TRANSFORMATIONS,
  const __global       int *SELECTIONS,
        __global PRECISION *TRANSFORMATIONS_BEST,
        __global       int *SELECTIONS_BEST,
  const               uint  NUM_COLS,
  const               uint  NUM_FRAGMENTS,
  const               uint  NUM_COORDINATES,
  const               uint  NUM_ROWS_MATRIX,
        __global       int *COUNTS_BEST,
        __global       int *DIFFERENCES,
        __local        int *SHARED_SUM_SELECTION,
        __local        int *SHARED_SUM_SELECTION_BEST
)
{
  const size_t num_cols_fragments     = get_global_size( 1);
  const size_t fragment_number_global = get_global_id( 1);
  const size_t fragment_number_local  = get_local_id( 1);
  const size_t number_fragments_local = get_local_size( 1);
  const size_t block_size             = get_local_size( 0);
  const size_t index_x                = get_local_id( 0);

  if( fragment_number_global >= NUM_FRAGMENTS)
  {
    return;
  }

  SHARED_SUM_SELECTION[ index_x * number_fragments_local + fragment_number_local] = 0;
  SHARED_SUM_SELECTION_BEST[ index_x * number_fragments_local + fragment_number_local] = 0;

  // reduce the selections
  for( size_t row = index_x; row < NUM_COORDINATES; row += block_size)
  {
    SHARED_SUM_SELECTION[ index_x * number_fragments_local + fragment_number_local]      += SELECTIONS[ row * num_cols_fragments + fragment_number_global];
    SHARED_SUM_SELECTION_BEST[ index_x * number_fragments_local + fragment_number_local] += SELECTIONS_BEST[ row * num_cols_fragments + fragment_number_global];
  }
  barrier( CLK_LOCAL_MEM_FENCE);

  // reduce accross threads
  for( size_t offset = block_size / 2; offset > 0; offset /= 2)
  {
    if( index_x < offset)
    {
      SHARED_SUM_SELECTION[ index_x * number_fragments_local + fragment_number_local]      += SHARED_SUM_SELECTION[ ( index_x + offset) * number_fragments_local + fragment_number_local];
      SHARED_SUM_SELECTION_BEST[ index_x * number_fragments_local + fragment_number_local] += SHARED_SUM_SELECTION_BEST[ ( index_x + offset) * number_fragments_local + fragment_number_local];
    }
    barrier( CLK_LOCAL_MEM_FENCE);
  }

  const int difference = SHARED_SUM_SELECTION[ fragment_number_local] - SHARED_SUM_SELECTION_BEST[ fragment_number_local];
  if( index_x == 0)
  {
    DIFFERENCES[ fragment_number_global] = difference;
  }

  // is best selection actually better than the best
  if( difference > 0)
  {
    if( index_x == 0)
    {
      COUNTS_BEST[ fragment_number_global] = SHARED_SUM_SELECTION[ fragment_number_local];
    }
    // update the best with the current selection
    for( size_t row = index_x; row < NUM_COORDINATES; row += block_size)
    {
      SELECTIONS_BEST[ row * num_cols_fragments + fragment_number_global] = SELECTIONS[ row * num_cols_fragments + fragment_number_global];
    }

    for( size_t copy_index = index_x; copy_index < NUM_ROWS_MATRIX * NUM_COLS; copy_index += block_size)
    {
      TRANSFORMATIONS_BEST[ fragment_number_global * NUM_ROWS_MATRIX * NUM_COLS + copy_index] = TRANSFORMATIONS[ fragment_number_global * NUM_ROWS_MATRIX * NUM_COLS + copy_index];
    }
  }
  else
  {
    if( index_x == 0)
    {
      COUNTS_BEST[ fragment_number_global] = SHARED_SUM_SELECTION_BEST[ fragment_number_local];
    }
  }
}
