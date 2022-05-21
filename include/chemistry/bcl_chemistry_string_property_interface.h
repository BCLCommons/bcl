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

#ifndef BCL_CHEMISTRY_STRING_PROPERTY_INTERFACE_H_
#define BCL_CHEMISTRY_STRING_PROPERTY_INTERFACE_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "util/bcl_util_enumerated.h"
#include "util/bcl_util_serializable_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class StringPropertyInterface
    //! @brief is a class to derive properties of small molecules or atoms
    //! @details it is derived from util::FunctionInterface to be able to use the framework around FunctionInterface
    //!
    //! @remarks example unnecessary
    //! @author mendenjl
    //! @date Oct 31, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    class BCL_API StringPropertyInterface :
      public util::SerializableInterface
    {
    public:

      //! virtual copy constructor
      virtual StringPropertyInterface *Clone() const = 0;

    ////////////////
    // operations //
    ////////////////

      //! @brief return a string given the small molecule
      //! @param MOLECULE the molecule to calculate the string property for
      //! @return the property as a string
      virtual std::string operator()( const ConformationInterface &MOLECULE) const = 0;

    }; // class StringPropertyInterface

  } // namespace chemistry
} // namespace bcl

#endif //BCL_CHEMISTRY_STRING_PROPERTY_INTERFACE_H_
