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

#ifndef BCL_QUALITY_DME_H_
#define BCL_QUALITY_DME_H_

// include the namespace header
#include "bcl_quality.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_quality_measure_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace quality
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class DME
    //! @brief Distance Matrix Error (DME) evaluates the difference between the distance matrices for given coordinates
    //! @details DME compares the distance matrices for two given coordinate sets of same size. This comparison is relatively
    //! fast since it does not require a superimposition of the coordinates beforehand.
    //!
    //! @see @link example_quality_dme.cpp @endlink
    //! @author karakam
    //! @date Oct 8, 2009
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API DME :
      public MeasureInterface
    {

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      DME();

      //! @brief Clone function
      //! @return pointer to new DME
      DME *Clone() const;

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
        return 0.0;
      }

      //! @brief return the comparison function for better quality
      //! @return binary function to compare two quality measure values
      const util::BinaryFunctionInterface< double, double, bool> &GetComparisonFunction() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief calculates distance matrix error (DME) between COORDINATES and REFERENCE_COORDINATES
      //! @param COORDINATES coordinates of interest
      //! @param REFERENCE_COORDINATES reference coordinates
      //! @return distance matrix error (DME) between COORDINATES and REFERENCE_COORDINATES
      double CalculateMeasure
      (
        const util::SiPtrVector< const linal::Vector3D> &COORDINATES,
        const util::SiPtrVector< const linal::Vector3D> &REFERENCE_COORDINATES
      ) const;

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

    }; // class DME

  } // namespace quality
} // namespace bcl

#endif // BCL_QUALITY_DME_H_ 
