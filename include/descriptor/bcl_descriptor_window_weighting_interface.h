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

#ifndef BCL_DESCRIPTOR_WINDOW_WEIGHTING_INTERFACE_H_
#define BCL_DESCRIPTOR_WINDOW_WEIGHTING_INTERFACE_H_

// include the namespace header
#include "bcl_descriptor.h"

// include other forward headers - sorted alphabetically
#include "biol/bcl_biol.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_descriptor_base_element.h"
#include "bcl_descriptor_window.h"
#include "util/bcl_util_function_interface_serializable.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class WindowWeightingInterface
    //! @brief given a half window size, generate a vector of weights
    //!
    //! @remarks example unnecessary
    //! @author mendenjl
    //! @date Mar 15, 2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API WindowWeightingInterface :
      public util::SerializableInterface
    {

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief Clone function
      //! @return pointer to new WindowWeightingInterface
      virtual WindowWeightingInterface *Clone() const = 0;

    /////////////////
    // data access //
    /////////////////

    ///////////////
    // operators //
    ///////////////

      //! @brief given a half window size, return the half window weights
      //! @param HALF_WINDOW_SIZE the size of the half window, counting the center element
      //! @return a vector of HALF_WINDOW_SIZE with the desired weights
      virtual linal::Vector< float> operator()( const size_t &HALF_WINDOW_SIZE) const = 0;

    }; // class WindowWeightingInterface

  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_WINDOW_WEIGHTING_INTERFACE_H_
