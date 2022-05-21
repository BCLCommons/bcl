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

#ifndef BCL_FIND_LOCATOR_COORDINATES_INTERFACE_H_
#define BCL_FIND_LOCATOR_COORDINATES_INTERFACE_H_

// include the namespace header
#include "bcl_find.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_find_locator_interface.h"
#include "linal/bcl_linal_vector_3d.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace find
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class LocatorCoordinatesInterface
    //! @brief Named interface for LocatorInterface
    //!
    //! @remarks example unnecessary
    //! @author mendenjl
    //! @date Nov 01, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    template< typename t_Argument>
    class LocatorCoordinatesInterface :
      public LocatorInterface< linal::Vector3D, t_Argument>
    {

    public:

      //! @brief Clone function
      //! @return pointer to new LocatorCoordinatesInterface
      virtual LocatorCoordinatesInterface *Clone() const = 0;

    }; // template class LocatorCoordinatesInterface

  } // namespace find
} // namespace bcl

#endif //BCL_FIND_LOCATOR_COORDINATES_INTERFACE_H_
