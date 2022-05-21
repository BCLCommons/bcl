// this is only included, so that eclipse c++ editor can be used and parses information correctly
// this will not be interpreted, when the program is compiled with an OpenCL compiler
#ifndef __OPENCL_VERSION__
 #include "kernel.h"
#endif

//! @brief Gathers k-th smallest values for each column of the matrix and inserts into the first k positions
//! @param INPUT        input matrix
//! @param INPUT_COLS   input cols including padding
//! @param INDEXES      index matrix
//! @param INDEXES_COLS index matrix cols including padding
//! @param COLS         total cols of input matrix
//! @param ROWS         total rows of input matrix
//! @param k            number of lowest values to find
__kernel void InsertionSort
(
  __global PRECISION* INPUT,
  __const  int        INPUT_COLS,
  __global int      * INDEXES,
  __const  int        INDEXES_COLS,
  __const  int        COLS,
  __const  int        ROWS,
  __const  int        k
)
{
  // Variables
  int l, i, j;
  __global PRECISION *p_dist;
  __global int   *p_ind;
  PRECISION curr_dist, max_dist;
  int   curr_row,  max_row;
  unsigned int xIndex = get_global_id( 0);

  if( xIndex < COLS)
  {
    // Pointer shift and initialization
    p_dist   = INPUT + xIndex;
    p_ind    = INDEXES  + xIndex;
    max_dist = p_dist[ 0];
    p_ind[ 0] = 0;

    // sort kth first element
    for ( l = 1; l < k; l++)
    {
      curr_row  = l * INPUT_COLS;
      curr_dist = p_dist[ curr_row];
      if( curr_dist < max_dist)
      {
        i = l - 1;
        for( int a = 0; a < l - 1; a++)
        {
          if( p_dist[ a * INPUT_COLS] > curr_dist)
          {
            i = a;
            break;
          }
        }
        for( j = l; j > i; j--)
        {
          p_dist[ j * INPUT_COLS] = p_dist[ ( j - 1) * INPUT_COLS];
          p_ind[ j * INDEXES_COLS]   = p_ind[ ( j - 1) * INDEXES_COLS];
        }
        p_dist[ i * INPUT_COLS] = curr_dist;
        p_ind[ i * INDEXES_COLS]   = l;// + 1;
      }
      else
      p_ind[ l * INDEXES_COLS] = l;// + 1;
      max_dist = p_dist[ curr_row];
    }

    // insert element in the k-th first lines
    max_row = ( k - 1) * INPUT_COLS;
    for( l = k; l < ROWS; l++)
    {
      curr_dist = p_dist[ l * INPUT_COLS];
      if( curr_dist < max_dist)
      {
        i = k - 1;
        for( int a = 0; a < k - 1; a++)
        {
          if( p_dist[ a * INPUT_COLS] > curr_dist)
          {
            i = a;
            break;
          }
        }
        for( j = k - 1; j > i; j--)
        {
          p_dist[ j * INPUT_COLS] = p_dist[ ( j - 1) * INPUT_COLS];
          p_ind[ j * INDEXES_COLS]   = p_ind[ ( j - 1) * INDEXES_COLS];
        }
        p_dist[ i * INPUT_COLS] = curr_dist;
        p_ind[ i * INDEXES_COLS]   = l;// + 1;
        max_dist                = p_dist[ max_row];
      }
    }
  }
}
