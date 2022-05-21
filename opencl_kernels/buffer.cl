// this is only included, so that eclipse c++ editor can be used and parses information correctly
// this will not be interpreted, when the program is compiled with an OpenCL compiler
#ifndef __OPENCL_VERSION__
 #include "kernel.h"
#endif

//! @brief fill all elements of buffer with a value
//! @param BUFFER the buffer to fill
//! @param OFFSET first index to fill
//! @param NUMBER_ELEMENTS the number of elements to fill
//! @param VALUE the value
__kernel void FillBuffer
(
        __global PRECISION* BUFFER,
  const               uint  OFFSET,
  const               uint  NUMBER_ELEMENTS,
  const          PRECISION  VALUE
)
{
  const size_t index = get_global_id( 0);
  if( index < NUMBER_ELEMENTS)
  {
    BUFFER[ index + OFFSET] = VALUE;
  }
}

//! @brief fill a submatrix
//! @param BUFFER the buffer to fill
//! @param VALUE the value
__kernel void FillSubmatrix
(
        __global PRECISION* BUFFER,
  const               uint  NUMBER_COLS,
  const               uint  ROW_OFFSET,
  const               uint  ROWS,
  const               uint  COL_OFFSET,
  const               uint  COLS,
  const          PRECISION  VALUE
)
{
  const size_t index_row = get_global_id( 0);
  const size_t index_col = get_global_id( 1);

  if( index_row < ROWS && index_col < COLS)
  {
    BUFFER[ ( index_row + ROW_OFFSET) * NUMBER_COLS + index_col + COL_OFFSET] = VALUE;
  }
}

//! @brief fill strided
//! @param BUFFER the buffer to fill
//! @param STRIDE offset between elements
//! @param NUMBER_ELEMENTS
//! @param VALUE the VALUE to assign to the cells
__kernel void FillStrided
(
        __global PRECISION* BUFFER,
  const               uint  STRIDE,
  const               uint  NUMBER_ELEMENTS,
  const          PRECISION  VALUE
)
{
  const size_t index = get_global_id( 0) * STRIDE;
  if( index < NUMBER_ELEMENTS)
  {
    BUFFER[ index] = VALUE;
  }
}

//! @brief fill an array with seeds like they are required by the GDT implementation
//! first row is filled starting from col 0 to seed length, second from 1 to seed_length +1 and so on
//! x dimension sets rows, y dimension each element within the seed at position x within row
//! @param BUFFER the buffer to fill
//! @param SEED_LENGTH length of seed
//! @param NUMBER_COLS number of cols in BUFFER
//! @param NUMBER_ROWS number of rows in BUFFER
//! @param VALUE the VALUE to assign to the cells
__kernel void FillSeeds
(
        __global PRECISION* BUFFER,
  const               uint  SEED_LENGTH,
  const               uint  NUMBER_COLS,
  const               uint  NUMBER_ROWS,
  const          PRECISION  VALUE
)
{
  const size_t seed_pos = get_global_id( 1);
  if( seed_pos >= SEED_LENGTH)
  {
    return;
  }
  const size_t row = get_global_id( 0);
  if( row + SEED_LENGTH > NUMBER_ROWS)
  {
    return;
  }

  // assign value
  BUFFER[ ( row + seed_pos) * NUMBER_COLS + row] = VALUE;
}
