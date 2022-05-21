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
#include "quality/bcl_quality_dmf.h"

// includes from bcl - sorted alphabetically
#include "coord/bcl_coord.h"
#include "linal/bcl_linal_matrix.h"
#include "linal/bcl_linal_vector_3d.h"
#include "math/bcl_math_comparisons.h"
#include "util/bcl_util_si_ptr_vector.h"

// external includes - sorted alphabetically
#include <numeric>

namespace bcl
{
  namespace quality
  {
  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> DMF::s_Instance
    (
      GetObjectInstances().AddInstance( new DMF( storage::Set< double>()))
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor from a distance cutoff set
    //! @param DISTANCE_CUTOFFS set of distance cutoffs
    DMF::DMF( const storage::Set< double> &DISTANCE_CUTOFFS) :
      m_DistanceCutoffs( DISTANCE_CUTOFFS)
    {
    }

    //! @brief Clone function
    //! @return pointer to new DMF
    DMF *DMF::Clone() const
    {
      return new DMF( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &DMF::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the comparison function for better quality
    //! @return binary function to compare two quality measure values
    const util::BinaryFunctionInterface< double, double, bool> &DMF::GetComparisonFunction() const
    {
      return **math::Comparisons< double>::GetEnums().e_Greater;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief calculates DMF between COORDINATES and REFERENCE_COORDINATES
    //! @param COORDINATES coordinates of interest
    //! @param REFERENCE_COORDINATES reference coordinates
    //! @return distance matrix fraction (DMF) between COORDINATES and REFERENCE_COORDINATES
    double DMF::CalculateMeasure
    (
      const util::SiPtrVector< const linal::Vector3D> &COORDINATES,
      const util::SiPtrVector< const linal::Vector3D> &REFERENCE_COORDINATES
    ) const
    {
      // make sure both coordinates are of same size
      BCL_Assert
      (
        COORDINATES.GetSize() == REFERENCE_COORDINATES.GetSize(),
        "The sizes of provided coordinates pair differ : " +
        util::Format()( COORDINATES.GetSize()) + " vs " + util::Format()( REFERENCE_COORDINATES.GetSize())
      );

      // calculate and store the difference of distance matrices
      const linal::Matrix< double> difference_matrix
      (
        coord::CalculateDifferenceDistanceMatrix( COORDINATES, REFERENCE_COORDINATES)
      );

      // dmf_sum
#define BCL_GCC_VERSION (__GNUC__ * 10000 + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__)
#if defined(__APPLE__) && BCL_GCC_VERSION < BCL_GCC_VERSION <= 40204
      double dmf_sum( 0.0);
      for
      (
        storage::Set< double>::const_iterator
          cutoff_itr( m_DistanceCutoffs.Begin()), cutoff_itr_end( m_DistanceCutoffs.End());
        cutoff_itr != cutoff_itr_end; ++cutoff_itr
      )
      {
        // sum up the value
        dmf_sum += CalculateMeasureMatrixCutoff( difference_matrix, *cutoff_itr);
      }
#else
      std::list< double> results_each_cutoff;
      std::transform
      (
        m_DistanceCutoffs.Begin(), m_DistanceCutoffs.End(),
        std::inserter( results_each_cutoff, results_each_cutoff.begin()),
        std::bind1st( std::ptr_fun( &CalculateMeasureMatrixCutoff), difference_matrix)
      );
      const double dmf_sum( std::accumulate( results_each_cutoff.begin(), results_each_cutoff.end(), double( 0)));
#endif

      // end
      return dmf_sum / double( m_DistanceCutoffs.GetSize());
    }

    //! @brief calculates DMF between COORDINATES and REFERENCE_COORDINATES for a given distance cutoff
    //! @param COORDINATES coordinates of interest
    //! @param REFERENCE_COORDINATES reference coordinates
    //! @param DISTANCE_CUTOFF distance cutoff
    //! @return DMF between COORDINATES and REFERENCE_COORDINATES for a given distance cutoff
    double DMF::CalculateMeasure
    (
      const util::SiPtrVector< const linal::Vector3D> &COORDINATES,
      const util::SiPtrVector< const linal::Vector3D> &REFERENCE_COORDINATES,
      const double DISTANCE_CUTOFF
    )
    {
      // make sure both coordinates are of same size
      BCL_Assert
      (
        COORDINATES.GetSize() == REFERENCE_COORDINATES.GetSize(),
        "The sizes of provided coordinates pair differ : " +
        util::Format()( COORDINATES.GetSize()) + " vs " + util::Format()( REFERENCE_COORDINATES.GetSize())
      );

      // calculate and store the difference of distance matrices
      const linal::Matrix< double> difference_matrix
      (
        coord::CalculateDifferenceDistanceMatrix( COORDINATES, REFERENCE_COORDINATES)
      );

      // calculate dmf and return it
      return CalculateMeasureMatrixCutoff( difference_matrix, DISTANCE_CUTOFF);
    }

    //! @brief calculates DMF for given difference matrix and a distance cutoff
    //! @param DIFFERENCE_MATRIX distance matrix that corresponds to two coordinate sets
    //! @param DISTANCE_CUTOFF  distance cutoff
    //! @return dmf for the given matrix
    double DMF::CalculateMeasureMatrixCutoff
    (
      const linal::Matrix< double> &DIFFERENCE_MATRIX,
      const double DISTANCE_CUTOFF
    )
    {
      // make sure the matrix is symmetrix
      BCL_Assert
      (
        DIFFERENCE_MATRIX.GetNumberRows() == DIFFERENCE_MATRIX.GetNumberCols(),
        "difference matrix is not symmetric: " + util::Format()( DIFFERENCE_MATRIX.GetNumberRows()) + " rows vs. " +
        util::Format()( DIFFERENCE_MATRIX.GetNumberCols()) + " cols!"
      );

      // initialize counts under the cutoff
      size_t counts( 0);

      // store number rows and columns
      const size_t number_rows( DIFFERENCE_MATRIX.GetNumberRows());

      // iterate over rows
      for( size_t i( 0); i < number_rows; ++i)
      {
        // iterate over the upper triangle
        for( size_t j( i + 1); j < number_rows; ++j)
        {
          // if the current cell has a value smaller than the cutoff
          if( DIFFERENCE_MATRIX( i, j) <= DISTANCE_CUTOFF)
          {
            // increment the count
            ++counts;
          }
        }
      }

      // do the normalization and return the value
      return double( 2.0) * DMF( storage::Set< double>()).OptimalValue() * double( counts) / double( number_rows * ( number_rows - 1));
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &DMF::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_DistanceCutoffs, ISTREAM);
      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &DMF::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_DistanceCutoffs, OSTREAM, INDENT);
      // end
      return OSTREAM;
    }

  } // namespace quality
} // namespace bcl
