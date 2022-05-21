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

#ifndef BCL_QUALITY_MEASURE_INTERFACE_H_
#define BCL_QUALITY_MEASURE_INTERFACE_H_

// include the namespace header
#include "bcl_quality.h"

// include other forward headers - sorted alphabetically
#include "linal/bcl_linal.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_object_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace quality
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MeasureInterface
    //! @brief MeasureInterface class provides the interface for different quality measures for comparing two coordinate sets
    //! @details MeasureInterface provides an operator that takes two coordinate vectors and returns a quality measure.
    //!
    //! @remarks example unnecessary
    //! @author alexanns, staritrd, woetzen, karakam
    //! @date 03/04/2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MeasureInterface :
      public util::ObjectInterface
    {

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief virtual copy constructor
      //! @return pointer to new Interface
      virtual MeasureInterface *Clone() const = 0;

    /////////////////
    // data access //
    /////////////////

      //! @brief return the optimal value for that quality measurement
      //! @return the best value by which two sets of coordinates can agree
      virtual double OptimalValue() const = 0;

      //! @brief return the comparison function for better quality
      //! @return binary function to compare two quality measure values
      virtual const util::BinaryFunctionInterface< double, double, bool> &GetComparisonFunction() const = 0;

    ////////////////
    // operations //
    ////////////////

      //! @brief calculates a quality measure between COORDINATES and REFERENCE_COORDINATES
      //! @param COORDINATES coordinates of interest
      //! @param REFERENCE_COORDINATES reference coordinates
      //! @return a double which is the agreement of COORDINATES between REFERENCE_COORDINATES
      virtual double CalculateMeasure
      (
        const util::SiPtrVector< const linal::Vector3D> &COORDINATES,
        const util::SiPtrVector< const linal::Vector3D> &REFERENCE_COORDINATES
      ) const = 0;

    }; // class MeasureInterface

  } // namespace quality
} // namespace bcl

#endif //BCL_QUALITY_MEASURE_INTERFACE_H_
