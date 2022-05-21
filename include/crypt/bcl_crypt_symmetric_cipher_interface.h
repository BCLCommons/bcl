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

#ifndef BCL_CRYPT_SYMMETRIC_CIPHER_INTERFACE_H_
#define BCL_CRYPT_SYMMETRIC_CIPHER_INTERFACE_H_

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
    //! @class SymmetricCipherInterface
    //! @brief interface to enforce common functionality among all symmetric cipher implementations
    //!
    //! @see no example necessary
    //! @author butkiem1
    //! @date Nov 3, 2011
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SymmetricCipherInterface :
      public util::ObjectInterface
    {

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

    /////////////////
    // data access //
    /////////////////

      //! @brief set password and related cryptographic tables
      //! @param PASSWORD password of interes
      virtual void SetPassword( const std::string &PASSWORD) = 0;

      //! @brief encrypt a plain text message and return an encrypted message
      //! @param ORIGINAL_MESSAGE plain text message as a string
      //! @return encrypted message
      virtual std::string Encrypt( const std::string &ORIGINAL_MESSAGE) = 0;

      //! @brief decrypt a encrypted message and return the original plain text
      //! @param ENCRYPTED_MESSAGE encrypted message as a string
      //! @return decrypted plain text
      virtual std::string Decrypt( const std::string &ENCRYPTED_MESSAGE) = 0;

    }; // class SymmetricCipherInterface

  } // namespace crypt
} // namespace bcl

#endif // BCL_CRYPT_SYMMETRIC_CIPHER_INTERFACE_H_ 
