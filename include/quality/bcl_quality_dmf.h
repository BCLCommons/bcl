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

#ifndef BCL_QUALITY_DMF_H_
#define BCL_QUALITY_DMF_H_

// include the namespace header
#include "bcl_quality.h"

// include other forward headers - sorted alphabetically
#include "linal/bcl_linal.fwd.hh"
#include "math/bcl_math.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_quality_measure_interface.h"
#include "storage/bcl_storage_set.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace quality
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class DMF
    //! @brief Distance Matrix Fraction evaluates the percentage fraction of elements in difference distance matrix that are under a given cutoff
    //! @details Distance Matrix Fraction (DMF) calculates the difference matrix between two coordinate sets of same size,
    //! similar to DME, but instead of calculating an average, it returns the ratio of elements in the matrix that
    //! are below a given distance cutoff. It also allows providing a vector of cutoffs and the return value is the
    //! average of the DMF values for each of the cutoffs.
    //!
    //! @see @link example_quality_dmf.cpp @endlink
    //! @author karakam
    //! @date Oct 8, 2009
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API DMF :
      public MeasureInterface
    {

    private:

    //////////
    // data //
    //////////

      //! vector of distance cutoffs
      storage::Set< double> m_DistanceCutoffs;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief constructor from a distance cutoff set
      //! @param DISTANCE_CUTOFFS set of distance cutoffs
      DMF( const storage::Set< double> &DISTANCE_CUTOFFS);

      //! @brief Clone function
      //! @return pointer to new DMF
      DMF *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief return the optimal value for that quality measurement
      //! @return the best value by which two sets of coordinates can agree
      double OptimalValue() const
      {
        return 100.0;
      }

      //! @brief return the comparison function for better quality
      //! @return binary function to compare two quality measure values
      const util::BinaryFunctionInterface< double, double, bool> &GetComparisonFunction() const;

      //! @brief get distance cutoffs
      //! @return distance cutoffs
      const storage::Set< double> &GetDistanceCutoffs() const
      {
        return m_DistanceCutoffs;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief calculates DMF between COORDINATES and REFERENCE_COORDINATES
      //! @param COORDINATES coordinates of interest
      //! @param REFERENCE_COORDINATES reference coordinates
      //! @return distance matrix fraction (DMF) between COORDINATES and REFERENCE_COORDINATES
      double CalculateMeasure
      (
        const util::SiPtrVector< const linal::Vector3D> &COORDINATES,
        const util::SiPtrVector< const linal::Vector3D> &REFERENCE_COORDINATES
      ) const;

      //! @brief calculates DMF between COORDINATES and REFERENCE_COORDINATES for a given distance cutoff
      //! @param COORDINATES coordinates of interest
      //! @param REFERENCE_COORDINATES reference coordinates
      //! @param DISTANCE_CUTOFF distance cutoff
      //! @return DMF between COORDINATES and REFERENCE_COORDINATES for a given distance cutoff
      static double CalculateMeasure
      (
        const util::SiPtrVector< const linal::Vector3D> &COORDINATES,
        const util::SiPtrVector< const linal::Vector3D> &REFERENCE_COORDINATES,
        const double DISTANCE_CUTOFF
      );

      //! @brief calculates DMF for given difference matrix and a distance cutoff
      //! @param DIFFERENCE_MATRIX distance matrix that corresponds to two coordinate sets
      //! @param DISTANCE_CUTOFF  distance cutoff
      //! @return dmf for the given matrix
      static double CalculateMeasureMatrixCutoff
      (
        const linal::Matrix< double> &DIFFERENCE_MATRIX,
        const double DISTANCE_CUTOFF
      );

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class DMF

  } // namespace quality
} // namespace bcl

#endif // BCL_QUALITY_DMF_H_
