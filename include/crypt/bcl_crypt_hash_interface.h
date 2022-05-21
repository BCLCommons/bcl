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

#ifndef BCL_CRYPT_HASH_INTERFACE_H_
#define BCL_CRYPT_HASH_INTERFACE_H_

// include the namespace header
#include "bcl_crypt.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "util/bcl_util_object_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace crypt
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class HashInterface
    //! @brief interface that can be used for hashing a string
    //!
    //! @see no example necessary
    //! @author mendenjl
    //! @date Mar 19, 2012
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API HashInterface :
      public util::ObjectInterface
    {

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

    /////////////////
    // data access //
    /////////////////

      //! @brief hash a plain text message and return an encrypted message
      //! @param ORIGINAL_MESSAGE plain text message as a string
      //! @return hashed message in format
      virtual std::string Hash( const std::string &ORIGINAL_MESSAGE) const = 0;

    }; // class HashInterface

  } // namespace crypt
} // namespace bcl

#endif // BCL_CRYPT_HASH_INTERFACE_H_ 
