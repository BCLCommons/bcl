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

#ifndef BCL_CRYPT_SHA1_H_
#define BCL_CRYPT_SHA1_H_

// include the namespace header
#include "bcl_crypt.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_crypt_hash_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace crypt
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Sha1
    //! @brief Secure hash standard - 1; computes a condensed representation of a message or a data file
    //! @details When a message of any length < 2^64 bits is input, the SHA-1 produces a 160-bit output called a message
    //!          digest. The message digest is very efficient to compute and makes comparison of very large strings
    //!          or files possible in near constant time because the message digest is usually much smaller in size than
    //!          the message.
    //!          @link http://www.itl.nist.gov/fipspubs/fip180-1.htm @endlink
    //!
    //! @see @link example_crypt_sha1.cpp @endlink
    //! @author mendenjl
    //! @date Mar 19, 2012
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Sha1 :
      public HashInterface
    {

    //////////
    // data //
    //////////

    private:

      // introducing enum with constants used throughout Sha1
      enum Constants
      {
        e_OutputSize = 20,                   //!< characters in the output
        e_OutputSizeInts = e_OutputSize / 4, //!< Ints in the output
        e_BlockSize = 64,                    //!< # of characters considered in each round
        e_WordBlockSize = 80,                //!< Size of working array
        e_MessagePadding = 8                 //!< Amount of padding reserved for size of string
      };

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief Clone function
      //! @return pointer to new Sha1
      Sha1 *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief encrypt a plain text message and return an encrypted message
      //! @param ORIGINAL_MESSAGE plain text message as a string
      //! @return hashed message (in hex)
      std::string Hash( const std::string &ORIGINAL_MESSAGE) const;

    ///////////////
    // operators //
    ///////////////

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

    //////////////////////
    // helper functions //
    //////////////////////

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //!
      //! @class Sha1State
      //! @brief Sha1State is a simple class that stores the state of the Sha1 algorithm
      //! @remarks example unnecessary
      //! @author mendenjl
      //! @date Mar 19, 2012
      //!
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      struct Sha1State
      {
        uint32_t m_MessageDigest[ e_OutputSizeInts]; //!< Message digest buffers

        //! @brief default constructor
        Sha1State( const std::string &STRING);

        //! @brief convert SHA-1 hashed string to hex string
        std::string ToHexString() const;

      private:

        //! @brief Process the next 512 bits of the message
        //! @param MESSAGE_BLOCK characters that make up the message
        //! @param STORAGE a temporary array used during processing
        void ProcessMessageBlock( unsigned char *MESSAGE_BLOCK, uint32_t *STORAGE);

      };

    }; // class Sha1

  } // namespace crypt
} // namespace bcl

#endif // BCL_CRYPT_SHA1_H_ 
