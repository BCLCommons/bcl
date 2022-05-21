// this is only included, so that eclipse c++ editor can be used and parses information correctly
// this will not be interpreted, when the program is compiled with an OpenCL compiler
#ifndef __OPENCL_VERSION__
 #include "kernel.h"
#endif

// 3 coordinates + 1 weight
#define NR_COORDINATE_AND_WEIGHT 4

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! @brief interpolates coordinates with their weight on a grid
//! @details the local wroksize is always 4x4x4, while the global worksize is the extension of the grid in x,y and z as
//!          multiple of 4
//! @param COORDINATES the coordinates to interpolate onto the grid, 3 elements for x, y and z plus the weight to apply
//! @param NR_COORDINATES nr of rows in COORDINATES
//! @param ORIGIN_X, the origin of the grid in real space x
//! @param ORIGIN_Y, the origin of the grid in real space y
//! @param ORIGIN_Z, the origin of the grid in real space z
//! @param WIDTH_X, the width of one grid element x
//! @param WIDTH_Y, the width of one grid element y
//! @param WIDTH_Z, the width of one grid element z
//! @param BLOB_K describes the shape of the gaussian blob
//! @param SD_CUTOFF square distance of atom to grid point, above which gaussian weight is not added
//! @param GRID the grid to which the coordinates and their weight are added interpolated to - grid spacing is 1.0
//! @param COORDINATES_SHARED matrix of shared coordinates
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
__kernel __attribute__((reqd_work_group_size(4, 4, 4))) void
SimulateDensityGaussianSphere
(
  __const __global PRECISION* COORDINATES,
  __const          uint       NR_COORDINATES,
  __const          PRECISION  ORIGIN_X,
  __const          PRECISION  ORIGIN_Y,
  __const          PRECISION  ORIGIN_Z,
  __const          PRECISION  WIDTH_X,
  __const          PRECISION  WIDTH_Y,
  __const          PRECISION  WIDTH_Z,
  __const          PRECISION  BLOB_K,
  __const          PRECISION  SD_CUTOFF,
          __global PRECISION* GRID,
          __local  PRECISION* COORDINATES_SHARED
)
{
  // global index
  __const size_t grid_x = get_global_id( 0);
  __const size_t grid_y = get_global_id( 1);
  __const size_t grid_z = get_global_id( 2);

  // work index
  __const size_t tx = get_local_id( 0);
  __const size_t ty = get_local_id( 1);
  __const size_t tz = get_local_id( 2);

  // intensity of current grid position
  PRECISION intensity = 0.0;

  // size of coordinate block copied to shared memory
  __const size_t sub_step = get_local_size( 0) * get_local_size( 1) * get_local_size( 2) * NR_COORDINATE_AND_WEIGHT;

  // reals space index of the grid element
  __const PRECISION real_space_x = ORIGIN_X + grid_x * WIDTH_X;
  __const PRECISION real_space_y = ORIGIN_Y + grid_y * WIDTH_Y;
  __const PRECISION real_space_z = ORIGIN_Z + grid_z * WIDTH_Z;

  // iterate over subsets of the coordinates to store them in shared memory
  for( size_t sub = 0; sub < NR_COORDINATES * NR_COORDINATE_AND_WEIGHT; sub += sub_step)
  {
    // copy 64 coordinates with 4 values each (3 corrected by the origin, the last is the weight - each work item copies one coordinate
    __const size_t index = tx * get_local_size( 0) * get_local_size( 1) * get_local_size( 2) + ty * get_local_size( 1) * get_local_size( 2) + tz * get_local_size( 2);

    // one work item copies all 4 elements of each coordinate
    COORDINATES_SHARED[ index + 0] = COORDINATES[ sub + index + 0];
    COORDINATES_SHARED[ index + 1] = COORDINATES[ sub + index + 1];
    COORDINATES_SHARED[ index + 2] = COORDINATES[ sub + index + 2];
    COORDINATES_SHARED[ index + 3] = COORDINATES[ sub + index + 3];

    // Synchronize to make sure that all coordinates are copied
    barrier( CLK_LOCAL_MEM_FENCE);

    // calculate distance of the sub_step coordinates to the grid position
    for( size_t coord = 0; coord < sub_step; coord += NR_COORDINATE_AND_WEIGHT)
    {
      __const PRECISION dx     = COORDINATES_SHARED[ coord + 0] - real_space_x;
      __const PRECISION dy     = COORDINATES_SHARED[ coord + 1] - real_space_y;
      __const PRECISION dz     = COORDINATES_SHARED[ coord + 2] - real_space_z;
      __const PRECISION weight = COORDINATES_SHARED[ coord + 3];
      __const PRECISION sd     = dx * dx + dy * dy + dz * dz;
      if( sd <= SD_CUTOFF)
      {
        intensity += weight * exp( -BLOB_K * sd);
      }
    }

    // Synchronize to make sure that work items are done before coordinates are replaced with new set
    barrier( CLK_LOCAL_MEM_FENCE);
  }

  // copy to global memory
  GRID[ grid_z * get_global_size( 0) * get_global_size( 1) + grid_y * get_global_size( 0) + grid_x] = intensity;
}
