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
#include "sdf/bcl_sdf_mdl_header.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_sdf_mdl_header.cpp
  //!
  //! @author mendenjl
  //! @date Feb 29, 2012
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleSdfMdlHeader :
    public ExampleInterface
  {
  public:

    ExampleSdfMdlHeader *Clone() const
    {
      return new ExampleSdfMdlHeader( *this);
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

      // check default constructor
      sdf::MdlHeader header_default;

      // check constructor from # atoms, # bonds
      sdf::MdlHeader header( 123, 456);

    /////////////////
    // data access //
    /////////////////

      // check get number atoms / bonds, and constructor from string
      BCL_ExampleCheck( header.IsValid(),       true);
      BCL_ExampleCheck( header.GetNumberAtoms(), 123);
      BCL_ExampleCheck( header.GetNumberBonds(), 456);

      // test get string
      BCL_ExampleCheck( header.ToMdlLine(), "123456  0  0  0  0  0  0  0  0999 V2000");

      // test set from string with a valid string
      const std::string test_string( "513876  0  0  0  0  0  0  0  0999 V2000");
      header.SetFromMdlLine( test_string, 0);
      BCL_ExampleIndirectCheck( header.ToMdlLine(), test_string, "SetFromString");
      BCL_ExampleIndirectCheck( header.IsValid(),   true,        "SetFromString");
      BCL_ExampleIndirectCheck( header.GetNumberAtoms(), 513, "SetFromString");
      BCL_ExampleIndirectCheck( header.GetNumberBonds(), 876, "SetFromString");

      const std::string invalid_test_string( "513876  0  0  0  0  0  0  0  0999 V3");
      header.SetFromMdlLine( invalid_test_string, 0);
      BCL_ExampleIndirectCheck( header.IsValid(), false, "SetFromString(" + invalid_test_string + ")");

      const std::string pre_v2000_test_string( "513876 ");
      header.SetFromMdlLine( pre_v2000_test_string, 0);
      BCL_ExampleIndirectCheck( header.IsValid(), false, "SetFromString(" + pre_v2000_test_string + ", 0)");
      header.SetFromMdlLine( pre_v2000_test_string, 3);
      BCL_ExampleIndirectCheck( header.IsValid(), true, "SetFromString(" + pre_v2000_test_string + ", 3)");

    //////////////////////
    // input and output //
    //////////////////////

      // test io
      header.SetFromMdlLine( test_string, 0); // have to set the header back to a valid mdl string first
      BCL_ExampleCheck( TestBCLObjectIOForSymmetry( header, sdf::MdlHeader()), true);
      BCL_ExampleCheck( TestBCLObjectOutputDiffers( header, sdf::MdlHeader()), true);

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleSdfMdlHeader

  const ExampleClass::EnumType ExampleSdfMdlHeader::s_Instance
  (
    GetExamples().AddEnum( ExampleSdfMdlHeader())
  );

} // namespace bcl
