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
#include "sdf/bcl_sdf_mdl_entry_types.h"
#include "sdf/bcl_sdf_mdl_line_types.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_sdf_mdl_entry_types.cpp
  //!
  //! @author butkiem1, alexanns
  //! @date May 14, 2011
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleSdfMdlEntryTypes :
    public ExampleInterface
  {
  public:

    ExampleSdfMdlEntryTypes *Clone() const
    {
      return new ExampleSdfMdlEntryTypes( *this);
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

    /////////////////
    // data access //
    /////////////////

      // check each of the default lines, which are created from the entry types
      BCL_ExampleCheck
      (
        sdf::GetDefaultLine( sdf::e_HeaderLine),
        "  0  0  0  0  0  0  0  0  0  0999 V2000"
      );
      BCL_ExampleCheck
      (
        sdf::GetDefaultLine( sdf::e_AtomLine),
        "    0.0000    0.0000    0.0000 X   0  0  0  0  0  0  0  0  0  0  0  0"
      );
      BCL_ExampleCheck( sdf::GetDefaultLine( sdf::e_BondLine), "  0  0  0  0  0  0  0");
      BCL_ExampleCheck( sdf::GetDefaultLine( sdf::e_TerminationLine), "$$$$");

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleSdfMdlEntryTypes

  const ExampleClass::EnumType ExampleSdfMdlEntryTypes::s_Instance
  (
    GetExamples().AddEnum( ExampleSdfMdlEntryTypes())
  );

} // namespace bcl
