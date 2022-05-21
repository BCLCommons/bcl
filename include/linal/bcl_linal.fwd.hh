// (c) Copyright BCL @ Vanderbilt University 2014
// (c) BCL Homepage: http://www.meilerlab.org/bclcommons
// (c) BCL Code Repository: https://github.com/BCLCommons/bcl
// (c)
// (c) The BioChemical Library (BCL) was originally developed by contributing members of the Meiler Lab @ Vanderbilt University.
// (c)
// (c) The BCL is now made available as an open-source software package distributed under the permissive MIT license,
// (c) developed and maintained by the Meiler Lab at Vanderbilt University and contributing members of the BCL Commons. 
// (c)
// (c) External code contributions to the BCL are welcome. Please visit the BCL Commons GitHub page for information on how you can contribute.
// (c) 
// (c) This file is part of the BCL software suite and is made available under the MIT license.
// (c)

#ifndef BCL_LINAL_FWD_HH_
#define BCL_LINAL_FWD_HH_

// include bcl_defines.h header
#include "bcl_defines.h"

// external includes - sorted alphabetically
#include <cstddef>

// This file contains forward declarations for the linal namespace
// This file is mostly automatically generated
// Developers may add typedefs and template default parameters
// all other changes will be removed the next time this file is generated
namespace bcl
{
  namespace linal
  {
  /////////////////////
  // regular classes //
  /////////////////////

    class DistanceGeometry;
    class Vector2D;
    class Vector3D;

  //////////////////////
  // template classes //
  //////////////////////

    template< typename t_DataType>
    class Matrix;

    template< typename t_DataType>
    class Matrix2x2;

    template< typename t_DataType>
    class Matrix3x3;

    template< typename t_DataType>
    class MatrixConstInterface;

    template< typename t_DataType>
    class MatrixConstReference;

    template< typename t_DataType>
    class MatrixInterface;

    template< typename t_DataType>
    class MatrixInversionCholesky;

    template< typename t_DataType>
    class MatrixInversionGaussJordan;

    template< typename t_DataType>
    class MatrixInversionInterface;

    template< typename t_DataType>
    class MatrixInversionMoorePenrose;

    template< typename t_DataType>
    class MatrixReference;

    template< typename t_DataType>
    class Operations;

    template< typename t_DataType>
    class OperationsCPU;

    template< typename t_DataType>
    class OperationsInterface;

    template< typename t_DataType>
    class PrincipalComponentAnalysis;

    template< typename t_DataType>
    class SymmetricEigenSolver;

    template< typename t_DataType>
    class Vector;

    template< typename t_DataType>
    class VectorConstInterface;

    template< typename t_DataType>
    class VectorConstReference;

    template< typename t_DataType>
    class VectorInterface;

    template< typename t_DataType, size_t t_N>
    class VectorND;

    template< typename t_DataType>
    class VectorReference;

  //////////////
  // typedefs //
  //////////////

  } // namespace linal
} // namespace bcl

#endif // BCL_LINAL_FWD_HH_
