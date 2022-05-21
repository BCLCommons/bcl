// this is only included, so that eclipse c++ editor can be used and parses information correctly
// this will not be interpreted, when the program is compiled with an OpenCL compiler
#ifndef __OPENCL_VERSION__
 #include "kernel.h"
#endif

#define BLOCK_SIZE 128
#define PRECALC_CONST_SQR 157.913670417
#define PRECALC_CONST 12.5663706144
#define SOLVENT_DENSITY 0.334

float ComputeFormFactor( float16 PARAMS, float Q_SQR, float H_FACTOR, float WATER_FACTOR)
{
  float ExcludedVolumeParameter       = PARAMS.s0;
  float HydrationShellParameter       = PARAMS.s1;
  float SolventAccessableSurfaceArea  = PARAMS.s2;
  float DisplacedSolventVolume        = PARAMS.s3;
  float BoundHydrogen                 = PARAMS.s5;
  float A1                            = PARAMS.s6;
  float A2                            = PARAMS.s7;
  float A3                            = PARAMS.s8;
  float A4                            = PARAMS.s9;
  float B1                            = PARAMS.sa;
  float B2                            = PARAMS.sb;
  float B3                            = PARAMS.sc;
  float B4                            = PARAMS.sd;
  float C                             = PARAMS.se;
  
  float form = 0;
  
  form += A1 * native_exp( -( B1 * Q_SQR / PRECALC_CONST_SQR));   
  form += A2 * native_exp( -( B2 * Q_SQR / PRECALC_CONST_SQR));   
  form += A3 * native_exp( -( B3 * Q_SQR / PRECALC_CONST_SQR));   
  form += A4 * native_exp( -( B4 * Q_SQR / PRECALC_CONST_SQR));   
  form += C;

  form += BoundHydrogen * H_FACTOR;
  
  return form - ExcludedVolumeParameter * DisplacedSolventVolume * SOLVENT_DENSITY * native_exp( Q_SQR * -native_powr( DisplacedSolventVolume, native_divide( 2.0f, 3.0f)) / PRECALC_CONST) + HydrationShellParameter * SolventAccessableSurfaceArea * WATER_FACTOR;
  
}

//! @brief calculate inner sums
//! @param RESIDUE_FORM_FACTORS form factors
//! @param DISTANCE_MATRIX pairwise distance matrix
//! @param Q array of q values           
//! @param INNER_SUM_MATRIX the output matrix
//! @param NUMBER_ATOMS number of atoms
//! @param NUMBER_Q number of q values
//! @param CURRENT_Q_INDEX current q being operated on
//! @param DISTANCES_MATRIX_PADDING padding for the distance matrix which is required for the pairwise kernel
__kernel void CalculateInnerSums
(
  const __global PRECISION* RESIDUE_FORM_FACTORS, 
  const __global PRECISION* DISTANCE_MATRIX,
  const __global PRECISION* Q,            
        __global PRECISION* INNER_SUM_MATRIX,
  const          uint       NUMBER_ATOMS,
  const          uint       NUMBER_Q,
  const          uint       CURRENT_Q_INDEX,
  const          uint       DISTANCES_MATRIX_PADDING
)
{
  if( get_global_id( 0) >= NUMBER_ATOMS) return;
  uint row = get_global_id( 0);
  
  uint col_start = row + 1;
  PRECISION residue_a = RESIDUE_FORM_FACTORS[ row * NUMBER_Q + CURRENT_Q_INDEX];
  PRECISION inner_sum = 0;
  PRECISION q = Q[ CURRENT_Q_INDEX];
    
  for( uint col = col_start; col < NUMBER_ATOMS; ++col)
  {
 	  PRECISION residue_b = RESIDUE_FORM_FACTORS[ col * NUMBER_Q + CURRENT_Q_INDEX];
 	
 	  PRECISION dist = DISTANCE_MATRIX[ row * ( NUMBER_ATOMS + DISTANCES_MATRIX_PADDING) + col];
 
    // float so that hardware level math functions can be used "native_*"
    // PRECISION x = q * dist;
    float x = q * dist;
    // fix x = 0 case
    // inner_sum += residue_b * ( sin( x)) / x;
    if( q == 0.0)
    {
      inner_sum += residue_b;
    }
    else
    {
      inner_sum += residue_b * native_sin( x) / x;
    }
  }
  
  INNER_SUM_MATRIX[ CURRENT_Q_INDEX * NUMBER_ATOMS + row] = 2.0 * residue_a * inner_sum + ( residue_a * residue_a);
}

__kernel void CalculateResFF
(
  __const __global float16 *PARAMS,
  __const float Q,
          __global float *OUTPUT,
  __const uint NR_ATOMS,
  __const float H_FACTOR,
  __const float WATER_FACTOR
)
{
  size_t i = get_global_id( 0);
  if( i >= NR_ATOMS) return;
  
  float q = Q;
  float q2 = q * q;
  float16 params = PARAMS[ i];
  float resff = ComputeFormFactor( params, q2, H_FACTOR, WATER_FACTOR);
    
  OUTPUT[ i] = resff;
}

//! @brief kernel for doing distances on the fly so we don't store it to reduce memory
//! @param COORDS float4 of atom coordinates with w=0
//! @param RES_FF residue form factors matrix of nr atoms x nr q-values
//! @param Q the q-values array
//! @param INNER_SUMS the partial sums matrix of size nr q-values x nr atoms - this is the output
//! @param NR_Q the number of q-values
//! @param NR_ATOMS the number of atoms
__kernel void InnerSumsOTFGPU
(
  __const __global float4 *COORDS,
  __const __global float  *RES_FF,
  __const          float   Q,
          __global float  *INNER_SUMS,
  __const          uint    NR_ATOMS
)
{
  __local float4 shared_coords[ BLOCK_SIZE];
  __local float shared_res_ff[ BLOCK_SIZE];
  size_t i = get_global_id( 0);

  if( i >= NR_ATOMS) return;

  float current_res_ff = RES_FF[ i];
  float x, y, z;


  float4 current_coords = COORDS[ i];
  float4 new_coords;
  size_t bid   = get_group_id( 0);
  size_t tid   = get_local_id( 0);
  size_t gsize = get_num_groups( 0);
  size_t bsize = get_local_size( 0);
  size_t j, k;
  float dist, tmp, sum = 0, res_b;
  float current_q = Q;

  for( j = 0; j < gsize; ++j)
  {
    if(  j * bsize + tid < NR_ATOMS)
    {
      shared_coords[ tid] = COORDS[ j * bsize + tid];
      shared_res_ff[ tid] = RES_FF[ j * bsize + tid];
    }
    barrier( CLK_LOCAL_MEM_FENCE);

    for( k = 0; k < bsize; ++k)
    {
      if( j * bsize + k >= NR_ATOMS || j * bsize + k <= i) continue;
      new_coords = shared_coords[ k];

      x = current_coords.x - new_coords.x;
      y = current_coords.y - new_coords.y;
      z = current_coords.z - new_coords.z;
      dist = native_sqrt( x*x + y*y + z*z);
      
	  
      res_b = shared_res_ff[ k];

      if(dist == 0)
      {
        sum += res_b;
        continue;
      }

      tmp = current_q * dist;
      sum += res_b * native_divide( native_sin( tmp), tmp);
    }
    barrier( CLK_LOCAL_MEM_FENCE);
  }

  INNER_SUMS[ i] = 2.0f * current_res_ff * sum + ( current_res_ff * current_res_ff);
}

//! @brief kernel for doing distances on the fly so we don't store it to reduce memory
//! @param COORDS float4 of atom coordinates with w=0
//! @param RES_FF residue form factors matrix of nr atoms x nr q-values
//! @param Q the q-values array
//! @param INNER_SUMS the partial sums matrix of size nr q-values x nr atoms - this is the output
//! @param NR_Q the number of q-values
//! @param NR_ATOMS the number of atoms
__kernel void InnerSumsOTFZeroGPU
(
  __const __global float4 *COORDS,
  __const __global float  *RES_FF,
  __const          float   Q,
          __global float  *INNER_SUMS,
  __const          uint    NR_ATOMS
)
{
  __local float shared_res_ff[ BLOCK_SIZE];
  size_t i = get_global_id( 0);

  if( i >= NR_ATOMS) return;

  float current_res_ff = RES_FF[ i];

  size_t bid   = get_group_id( 0);
  size_t tid   = get_local_id( 0);
  size_t gsize = get_num_groups( 0);
  size_t bsize = get_local_size( 0);
  size_t j, k;
  float sum = 0, res_b;

  for( j = 0; j < gsize; ++j)
  {
    if(  j * bsize + tid < NR_ATOMS)
    {
      shared_res_ff[ tid] = RES_FF[ j * bsize + tid];
    }
    barrier( CLK_LOCAL_MEM_FENCE);

    for( k = 0; k < bsize; ++k)
    {
      if( j * bsize + k >= NR_ATOMS || j * bsize + k <= i) continue;

      res_b = shared_res_ff[ k];

      sum += res_b;
    }
    barrier( CLK_LOCAL_MEM_FENCE);
  }

  INNER_SUMS[ i] = 2.0f * current_res_ff * sum + ( current_res_ff * current_res_ff);
}

//! @brief kernel for doing distances on the fly so we don't store it to reduce memory
//! @param COORDS float4 of atom coordinates with w=0
//! @param RES_FF residue form factors matrix of nr atoms x nr q-values
//! @param Q the q-values array
//! @param INNER_SUMS the partial sums matrix of size nr q-values x nr atoms - this is the output
//! @param NR_Q the number of q-values
//! @param NR_ATOMS the number of atoms
__kernel void InnerSumsOTFCPU
(
  __const __global float4 *COORDS,
  __const __global float  *RES_FF,
  __const          float   Q,
          __global float  *INNER_SUMS,
  __const          uint    NR_ATOMS
)
{
  size_t i = get_global_id( 0);

  if( i >= NR_ATOMS) return;

  float current_res_ff = RES_FF[ i];
  float x, y, z;


  float4 current_coords = COORDS[ i];
  float4 new_coords;
  size_t bid   = get_group_id( 0);
  size_t tid   = get_local_id( 0);
  size_t gsize = get_num_groups( 0);
  size_t bsize = get_local_size( 0);
  size_t j, k;
  float dist, tmp, sum = 0, res_b;
  float current_q = Q;

  for( j = 0; j < gsize; ++j)
  {
    for( k = 0; k < bsize; ++k)
    {
      if( j * bsize + k >= NR_ATOMS || j * bsize + k <= i) continue;
      new_coords = COORDS[ j * bsize + k];
      res_b =  RES_FF[ j * bsize + k];

      x = current_coords.x - new_coords.x;
      y = current_coords.y - new_coords.y;
      z = current_coords.z - new_coords.z;

      dist = native_sqrt( x*x + y*y + z*z);

      tmp = current_q * dist;
      sum += res_b * native_divide( native_sin( tmp), tmp);
    }
  }

  INNER_SUMS[ i] = 2.0f * current_res_ff * sum + ( current_res_ff * current_res_ff);
}

//! @brief kernel for doing distances on the fly so we don't store it to reduce memory
//! @param COORDS float4 of atom coordinates with w=0
//! @param RES_FF residue form factors matrix of nr atoms x nr q-values
//! @param Q the q-values array
//! @param INNER_SUMS the partial sums matrix of size nr q-values x nr atoms - this is the output
//! @param NR_Q the number of q-values
//! @param NR_ATOMS the number of atoms
__kernel void InnerSumsOTFZeroCPU
(
  __const __global float4 *COORDS,
  __const __global float  *RES_FF,
  __const          float   Q,
          __global float  *INNER_SUMS,
  __const          uint    NR_ATOMS
)
{
  size_t i = get_global_id( 0);

  if( i >= NR_ATOMS) return;

  float current_res_ff = RES_FF[ i];

  size_t bid   = get_group_id( 0);
  size_t tid   = get_local_id( 0);
  size_t gsize = get_num_groups( 0);
  size_t bsize = get_local_size( 0);
  size_t j, k;
  float sum = 0, res_b;

  for( j = 0; j < gsize; ++j)
  {
    for( k = 0; k < bsize; ++k)
    {
      if( j * bsize + k >= NR_ATOMS || j * bsize + k <= i) continue;

      res_b =  RES_FF[ j * bsize + k];

      sum += res_b;
    }
  }

  INNER_SUMS[ i] = 2.0f * current_res_ff * sum + ( current_res_ff * current_res_ff);
}

//! @brief performs a column reduction
//! @param INPUT the input matrix of intensities
//! @param NUM_ELEMENTS number of elements
//! @param
//! @param
//! @param
__kernel void ReductionSumOfIntensity
(
  __const __global PRECISION* INPUT,
  __const          uint       NUM_ELEMENTS,
          __global PRECISION* OUTPUT,
          __local  PRECISION* SHARED
)
{
  // perform first level of reduction,
  // reading from global memory, writing to shared memory
  __const size_t thread_id = get_local_id( 0);
  __const size_t block_size = get_local_size( 0);
  __const size_t grid_size = block_size * 2 * get_num_groups( 0);
  size_t index = get_group_id( 0) * ( get_local_size( 0) * 2) + get_local_id( 0);
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
  if( block_size >= 512)
  {
    if( thread_id < 256)
    {
      SHARED[ thread_id] += SHARED[ thread_id + 256];
    }
    barrier( CLK_LOCAL_MEM_FENCE);
  }
  if( block_size >= 256)
  {
    if( thread_id < 128)
    {
      SHARED[ thread_id] += SHARED[ thread_id + 128];
    }
    barrier( CLK_LOCAL_MEM_FENCE);
  }
  if( block_size >= 128)
  {
    if( thread_id <  64)
    {
      SHARED[ thread_id] += SHARED[ thread_id +  64];
    }
    barrier( CLK_LOCAL_MEM_FENCE);
  }
  if( block_size >= 64)
  {
    if( thread_id <  32)
    {
      SHARED[ thread_id] += SHARED[ thread_id +  32];
    }
    barrier( CLK_LOCAL_MEM_FENCE);
  }
  if( block_size >= 32)
  {
    if( thread_id <  16)
    {
      SHARED[ thread_id] += SHARED[ thread_id +  16];
    }
    barrier( CLK_LOCAL_MEM_FENCE);
  }
  if( block_size >= 16)
  {
    if( thread_id <  8)
    {
      SHARED[ thread_id] += SHARED[ thread_id +  8];
    }
    barrier( CLK_LOCAL_MEM_FENCE);
  }
  if( block_size >= 8)
  {
    if( thread_id <  4)
    {
      SHARED[ thread_id] += SHARED[ thread_id +  4];
    }
    barrier( CLK_LOCAL_MEM_FENCE);
  }
  if( block_size >= 4)
  {
    if( thread_id <  2)
    {
      SHARED[ thread_id] += SHARED[ thread_id +  2];
    }
    barrier( CLK_LOCAL_MEM_FENCE);
  }
  if( block_size >= 2)
  {
    if( thread_id <  1)
    {
      SHARED[ thread_id] += SHARED[ thread_id +  1];
    }
    barrier( CLK_LOCAL_MEM_FENCE);
  }

  // write result for this block to global mem
  if( thread_id == 0)
  {
    OUTPUT[ get_group_id( 0)] = SHARED[ 0];
  }
  barrier( CLK_LOCAL_MEM_FENCE);
}
