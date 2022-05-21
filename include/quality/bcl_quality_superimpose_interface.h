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

#ifndef BCL_QUALITY_SUPERIMPOSE_INTERFACE_H_
#define BCL_QUALITY_SUPERIMPOSE_INTERFACE_H_

// include the namespace header
#include "bcl_quality.h"

// include other forward headers - sorted alphabetically
#include "math/bcl_math.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_quality_measure_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace quality
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SuperimposeInterface
    //! @brief Interface class for calculating a superimposition between two coordinate sets
    //! @details SuperimposeInterface provides a virtual Superimpose function that calculates the
    //! transformation matrix that gives the best superimposition between two coordinates according to a quality measure
    //!
    //! @remarks example unnecessary
    //! @author karakam
    //! @date Oct 31, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SuperimposeInterface :
      public MeasureInterface
    {

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief Clone function
      //! @return pointer to new SuperimposeInterface
      virtual SuperimposeInterface *Clone() const = 0;

    ////////////////
    // operations //
    ////////////////

      //! @brief calculates the transformation matrix that superimposes COORDINATES onto REFERENCE_COORDINATES
      //! @param COORDINATES Vector of coordinates of interest
      //! @param REFERENCE_COORDINATES reference coordinates
      //! @return transformation matrix that superimposes COORDINATES to REFERENCE_COORDINATES
      virtual math::TransformationMatrix3D CalculateSuperimposition
      (
        const util::SiPtrVector< const linal::Vector3D> &COORDINATES,
        const util::SiPtrVector< const linal::Vector3D> &REFERENCE_COORDINATES
      ) const = 0;

    }; // class SuperimposeInterface

  } // namespace quality
} // namespace bcl

#endif // BCL_QUALITY_SUPERIMPOSE_INTERFACE_H_ 
