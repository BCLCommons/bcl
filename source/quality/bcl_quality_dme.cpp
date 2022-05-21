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
#include "quality/bcl_quality_dme.h"

// includes from bcl - sorted alphabetically
#include "coord/bcl_coord.h"
#include "linal/bcl_linal_matrix_operations.h"
#include "linal/bcl_linal_vector_3d.h"
#include "math/bcl_math_comparisons.h"
#include "util/bcl_util_si_ptr_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace quality
  {
  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> DME::s_Instance
    (
      GetObjectInstances().AddInstance( new DME())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    DME::DME()
    {
    }

    //! @brief Clone function
    //! @return pointer to new DME
    DME *DME::Clone() const
    {
      return new DME( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &DME::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the comparison function for better quality
    //! @return binary function to compare two quality measure values
    const util::BinaryFunctionInterface< double, double, bool> &DME::GetComparisonFunction() const
    {
      return **math::Comparisons< double>::GetEnums().e_Less;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief calculates distance matrix error (DME) between COORDINATES and REFERENCE_COORDINATES
    //! @param COORDINATES coordinates of interest
    //! @param REFERENCE_COORDINATES reference coordinates
    //! @return distance matrix error (DME) between COORDINATES and REFERENCE_COORDINATES
    double DME::CalculateMeasure
    (
      const util::SiPtrVector< const linal::Vector3D> &COORDINATES,
      const util::SiPtrVector< const linal::Vector3D> &REFERENCE_COORDINATES
    ) const
    {
      BCL_Assert
      (
        COORDINATES.GetSize() == REFERENCE_COORDINATES.GetSize(),
        "number of given coordinates are different: " +
        util::Format()( COORDINATES.GetSize()) + " != " + util::Format()( REFERENCE_COORDINATES.GetSize())
      );

      // calculate and store the difference of distance matrices
      const linal::Matrix< double> difference_matrix
      (
        coord::CalculateDifferenceDistanceMatrix( COORDINATES, REFERENCE_COORDINATES)
      );

      // initialize dme value to calculate
      double dme_squared( 0.0);

      // make sure the matrix is symmetrix
      BCL_Assert
      (
        difference_matrix.IsSquare(), "difference matrix is not square matrix"
      );

      // store number rows and columns
      const size_t number_rows( difference_matrix.GetNumberRows());

      // iterate over rows
      for( size_t i( 0); i < number_rows; ++i)
      {
        // iterate over the upper triangle
        for( size_t j( i + 1); j < number_rows; ++j)
        {
          // sum up the square of the distance
          dme_squared += math::Sqr( difference_matrix( i, j));
        }
      }
      // do the normalization
      dme_squared *= 2.0 / double( number_rows * ( number_rows - 1));

      // return square root of the calculated value
      return math::Sqrt( dme_squared);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &DME::Read( std::istream &ISTREAM)
    {
      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &DME::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // end
      return OSTREAM;
    }

  } // namespace quality
} // namespace bcl
