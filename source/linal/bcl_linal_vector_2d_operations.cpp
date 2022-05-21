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

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "linal/bcl_linal_vector_2d_operations.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace linal
  {
    //! @brief calculates the absolute difference between two linal::Vector2Ds
    //! @param VECTOR_A first vector
    //! @param VECTOR_B second vector
    //! @return absolute difference between two linal::Vector2D (points)
    Vector2D AbsoluteDifference( const Vector2D &VECTOR_A, const Vector2D &VECTOR_B)
    {
      Vector2D returned_vector( VECTOR_B - VECTOR_A);
      math::Absolute( returned_vector);
      return returned_vector;
    }

    //! @brief calculates the unit vector starting from one linal::Vector2D to another
    //! @param ORIGIN vector of origin
    //! @param TARGET target vector
    //! @return the unit vector between ORIGIN and TARGET
    Vector2D UnitVector( const Vector2D &ORIGIN, const Vector2D &TARGET)
    {
      Vector2D unit_vector( TARGET);
      unit_vector -= ORIGIN;
      unit_vector.Normalize();
      return unit_vector;
    }

  } // namespace linal
} // namespace bcl
