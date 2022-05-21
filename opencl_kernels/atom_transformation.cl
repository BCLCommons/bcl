// this is only included, so that eclipse c++ editor can be used and parses information correctly
// this will not be interpreted, when the program is compiled with an OpenCL compiler
#ifndef __OPENCL_VERSION__
 #include "kernel.h"
#endif

// block size is 4 * 4
#define NR_COORDINATE_AND_PROPERTY 4 // 3 coordinates + 1 property
#define AS(i, j) ATOM_SHARED[ i * NR_COORDINATE_AND_PROPERTY + j]
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! @brief coordinate transformation of atoms on the device, one atoms has 3 coordinates and one additional property,
//!        that is ignored
//! @details it basically is a matrix multiplication NEW_ATOMS = ATOMS * TRANSFORMATION, only that the last col of atoms is treated as 1
//! @param NEW_ATOMS the new atoms as many as the ATOMS
//! @param ATOMS the current atoms with NR_ATOMS rows and NR_COORDINATE_AND_PROPERTY cols
//! @param TRANSFORMATION_MATRIX a 4x4  matrix that contains the translation in
//!        the last row, and 0 in the last col, except for the last element, which is 1
//! @param ATOM_SHARED shared memory
//! @param TRANSFORMATION_SHARED shared memory
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
__kernel void
CoordinateTransformation
(
  __const __global PRECISION* ATOMS,
  __const          PRECISION  A1,
  __const          PRECISION  A2,
  __const          PRECISION  A3,
  __const          PRECISION  B1,
  __const          PRECISION  B2,
  __const          PRECISION  B3,
  __const          PRECISION  C1,
  __const          PRECISION  C2,
  __const          PRECISION  C3,
  __const          PRECISION  T1,
  __const          PRECISION  T2,
  __const          PRECISION  T3,
          __global PRECISION* NEW_ATOMS
)
{
  // Thread index
  __const size_t atom = get_global_id( 0);

  // old x y and z coordinate
  __const PRECISION old_x = ATOMS[ atom * NR_COORDINATE_AND_PROPERTY + 0];
  __const PRECISION old_y = ATOMS[ atom * NR_COORDINATE_AND_PROPERTY + 1];
  __const PRECISION old_z = ATOMS[ atom * NR_COORDINATE_AND_PROPERTY + 2];

  // new x y and z coordinate
  PRECISION new_x = 0.0;
  PRECISION new_y = 0.0;
  PRECISION new_z = 0.0;

  // Multiply the two matrices together; each thread computes one element
  // of the block sub-matrix
  new_x += old_x * A1;
  new_x += old_y * B1;
  new_x += old_z * C1;
  new_x += T1;
  new_y += old_x * A2;
  new_y += old_y * B2;
  new_y += old_z * C2;
  new_y += T2;
  new_z += old_x * A3;
  new_z += old_y * B3;
  new_z += old_z * C3;
  new_z += T3;

  // Write the block sub-matrix to device memory;
  // each thread writes one element
  NEW_ATOMS[ atom * NR_COORDINATE_AND_PROPERTY + 0] = new_x;
  NEW_ATOMS[ atom * NR_COORDINATE_AND_PROPERTY + 1] = new_y;
  NEW_ATOMS[ atom * NR_COORDINATE_AND_PROPERTY + 2] = new_z;
  NEW_ATOMS[ atom * NR_COORDINATE_AND_PROPERTY + 3] = ATOMS[ atom * NR_COORDINATE_AND_PROPERTY + 3];
}
