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
#include "chemistry/bcl_chemistry_solvents.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_chemistry_solvents.cpp
  //!
  //! @author riddeljs
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleChemistrySolvents :
    public ExampleInterface
  {
  public:

    ExampleChemistrySolvents *Clone() const
    {
      return new ExampleChemistrySolvents( *this);
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

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

      // check IsDefined function, see if undefined solvent is indeed undefined
      BCL_ExampleCheck( chemistry::IsDefined( chemistry::e_UNDEFINED_SOLVENT), false);

      // check IsDefined function, see if a solven in the enum is actually defined
      BCL_ExampleCheck( chemistry::IsDefined( chemistry::e_DMSO), true);

      // test SolventFromString, looking up METHANOL should return METHANOL_D4
      BCL_ExampleCheck( chemistry::SolventFromString( "METHANOL"), chemistry::e_METHANOL_D4);

      // test SolventFromString, this one should return UNDEFINED_SOLVENT
      BCL_ExampleCheck( chemistry::SolventFromString( "METHANOLS"), chemistry::e_UNDEFINED_SOLVENT);

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleChemistrySolvents

  const ExampleClass::EnumType ExampleChemistrySolvents::s_Instance
  (
    GetExamples().AddEnum( ExampleChemistrySolvents())
  );

} // namespace bcl
