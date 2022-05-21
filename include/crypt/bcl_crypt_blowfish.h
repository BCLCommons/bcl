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

#ifndef BCL_CRYPT_BLOWFISH_H_
#define BCL_CRYPT_BLOWFISH_H_

// include the namespace header
#include "bcl_crypt.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_crypt_symmetric_cipher_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace crypt
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class BlowFish
    //! @brief BlowFish is a symmetric cipher that uses a password to encrypt and decrypt messages
    //! @details Blowfish is a keyed, symmetric block cipher, designed in 1993 by Bruce Schneier and included in a large
    //!          number of cipher suites and encryption products. Blowfish provides a good encryption rate in software
    //!          and no effective cryptanalysis of it has been found to date.
    //!          Blowfish has a 64-bit block size and a variable key length from 1 bit up to 448 bits.
    //!          @link http://en.wikipedia.org/wiki/Blowfish_%28cipher%29 @endlink
    //!
    //! @see @link example_crypt_blowfish.cpp @endlink
    //! @author butkiem1
    //! @date Nov 3, 2011
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API BlowFish :
      public SymmetricCipherInterface
    {

    //////////
    // data //
    //////////

    private:

      // introducing enum with constants used throughout blowfish
      enum Constants
      {
        number_subkeys = 18,
        number_substitution_boxes = 4,
        number_entries = 256,
        max_password_length = 56
      };

      // contains padding
      unsigned int PA[ number_subkeys];

      // contains substitution boxes of symmetric key algorithm
      unsigned int SB[ number_substitution_boxes][ number_entries];

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      BlowFish();

      //! @brief constructor with parameter
      BlowFish( const std::string &PASSWORD);

      //! @brief Clone function
      //! @return pointer to new BlowFish
      BlowFish *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief set password and related cryptographic tables
      //! @param PASSWORD password of interes
      void SetPassword( const std::string &PASSWORD);

      //! @brief encrypt a plain text message and return an encrypted message
      //! @param ORIGINAL_MESSAGE plain text message as a string
      //! @return encrypted message
      std::string Encrypt( const std::string &ORIGINAL_MESSAGE);

      //! @brief decrypt a encrypted message and return the original plain text
      //! @param ENCRYPTED_MESSAGE encrypted message as a string
      //! @return decrypted plain text
      std::string Decrypt( const std::string &ENCRYPTED_MESSAGE);

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

    private:

      //! @brief reset the cryptographic tables of member variables
      void Reset();

      //! @brief generate sub keys stored to update cryptographic tables based on given password
      //! @param PASSWORD given password as string
      void GenerateSubkeys( const std::string &PASSWORD);

      //! @brief encrypt according to the blowfish encryption algorithm
      //! @param BYTES_HIGH first four bytes of the message being encrypted
      //! @param BYTES_LOW next four bytes of the message being encrypted
      void BlowFishEncrypt( unsigned int &BYTES_HIGH, unsigned int &BYTES_LOW);

      //! @brief decrypt according to the blowfish encryption algorithm
      //! @param BYTES_HIGH first four bytes of the message being decrypted
      //! @param BYTES_LOW next four bytes of the message being decrypted
      void BlowFishDecrypt( unsigned int &BYTES_HIGH, unsigned int &BYTES_LOW);

      //! @brief core cipher function used in the blowfish encryption
      //! @param VALUE value the cipher is applied on
      unsigned int Cipher( const unsigned int VALUE);

      //! @brief convenience function to write a 1byte portion of a 4byte integer at a specific byte position
      //! @param SET integer of interest
      //! @param VALUE character beeing written onto integer
      //! @param POSITION byte position in integer
      void Assign( unsigned int &SET, const char &VALUE, const size_t POSITION) const;

      //! @brief convenience function to get a 1byte portion of a 4byte integer at a specific byte position
      //! @param SET integer of interest
      //! @param POSITION byte position in integer
      unsigned int Get( const unsigned int &SET, const size_t POSITION) const;

    }; // class BlowFish

  } // namespace crypt
} // namespace bcl

#endif // BCL_CRYPT_BLOWFISH_H_ 
