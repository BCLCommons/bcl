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
#include "crypt/bcl_crypt_sha1.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_crypt_sha1.cpp
  //! @brief example for symmetric cypher implementation Sha1.
  //!
  //! @author mendenjl
  //! @date Mar 19, 2012
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleCryptSha1 :
    public ExampleInterface
  {
  public:

    ExampleCryptSha1 *Clone() const
    {
      return new ExampleCryptSha1( *this);
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
      crypt::Sha1 cypher_default;

      // test clone
      util::ShPtr< crypt::Sha1> sp_cypher( cypher_default.Clone());

    /////////////////
    // data access //
    /////////////////

    ////////////////
    // operations //
    ////////////////

      // perform the examples checks listed in FIP 180-1
      //! @see http://www.itl.nist.gov/fipspubs/fip180-1.htm
      BCL_ExampleCheck( cypher_default.Hash( "abc"), "A9993E364706816ABA3E25717850C26C9CD0D89D");
      BCL_ExampleCheck
      (
        cypher_default.Hash( "abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq"),
        "84983E441C3BD26EBAAE4AA1F95129E5E54670F1"
      );
      BCL_ExampleCheck
      (
        cypher_default.Hash( std::string( 1000000, 'a')),
        "34AA973CD4C4DAA4F61EEB2BDBAD27316534016F"
      );

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

  }; //end ExampleCryptSha1

  const ExampleClass::EnumType ExampleCryptSha1::s_Instance
  (
    GetExamples().AddEnum( ExampleCryptSha1())
  );

} // namespace bcl
