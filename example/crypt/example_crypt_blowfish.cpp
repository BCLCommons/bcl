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

// include example header
#include "example.h"
// include the header of the class which this example is for
#include "crypt/bcl_crypt_blowfish.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_crypt_blowfish.cpp
  //! @brief example for symmetric cypher implementation Blowfish.
  //!
  //! @author butkiem1
  //! @date Dec 10, 2011
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleCryptBlowfish :
    public ExampleInterface
  {
  public:

    ExampleCryptBlowfish *Clone() const
    {
      return new ExampleCryptBlowfish( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    int Run() const
    {
    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      crypt::BlowFish cypher_default;

      // constructor with parameter
      crypt::BlowFish cypher( "pass");

      // test clone
      util::ShPtr< crypt::BlowFish> sp_cypher( cypher.Clone());

    /////////////////
    // data access //
    /////////////////

      // set a new password
      cypher.SetPassword( "other_password");

    ////////////////
    // operations //
    ////////////////

      // define a message string
      std::string message( "my_message");

      std::string encrypted_message, decrypted_message;

      // encrypt message string, note: encrypted message must have a length with a multiple of eight
      // if not the message gets adjusted in the encryption method!
      encrypted_message = cypher.Encrypt( message);

      BCL_ExampleCheck( encrypted_message.length(), size_t( 16));

      // decrypt encrypted string
      decrypted_message = cypher.Decrypt( encrypted_message);

      // remove space characters
      decrypted_message = util::RemoveSpacesFromString( decrypted_message);

      // check for length of decrypted string
      BCL_ExampleCheck( decrypted_message.length(), size_t( 10));

      // check the decrypted message
      BCL_ExampleCheck( decrypted_message == "my_message", true);

    ///////////////
    // operators //
    ///////////////

    //////////////////////
    // input and output //
    //////////////////////

    ////////////////////
    // read and write //
    ////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleCryptBlowfish

  const ExampleClass::EnumType ExampleCryptBlowfish::s_Instance
  (
    GetExamples().AddEnum( ExampleCryptBlowfish())
  );

} // namespace bcl
