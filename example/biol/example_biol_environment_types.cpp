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
#include "biol/bcl_biol_environment_types.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_biol_environment_types.cpp
  //!
  //! @author karakam
  //! @date
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on 
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleBiolEnvironmentTypes :
    public ExampleInterface
  {
  public:

    ExampleBiolEnvironmentTypes *Clone() const
    {
      return new ExampleBiolEnvironmentTypes( *this);
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

      // no constructor for EnvironemntTypes available

      // construct EnvironmentType
      biol::EnvironmentType default_env_type;

      // construct from name
      biol::EnvironmentType env_type_solution( "SOLUTION");

      // from entype
      biol::EnvironmentType env_type_membrane( biol::GetEnvironmentTypes().e_MembraneCore);

    /////////////////
    // data access //
    /////////////////

      // iterator
      biol::EnvironmentTypes::const_iterator itr_beg( biol::GetEnvironmentTypes().Begin());
      biol::EnvironmentTypes::const_iterator itr_end( biol::GetEnvironmentTypes().End());

      // check if there are environment types available
      BCL_ExampleCheck( itr_beg < itr_end, true);
      BCL_ExampleCheck( biol::GetEnvironmentTypes().GetEnumCount() != 0, true);

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

      // access to reduced types
      BCL_ExampleCheck( biol::GetEnvironmentTypes().GetReducedTypes().GetSize(), 3);

      // create type from two letter code
      const biol::EnvironmentType two_letter_code( biol::GetEnvironmentTypes().EnvironmentTypeFromTwoLetterCode( "TR"));
      BCL_ExampleCheck( two_letter_code, biol::GetEnvironmentTypes().e_Transition);

    //////////////////////
    // input and output //
    //////////////////////

      // write all enums
      WriteBCLObject( biol::GetEnvironmentTypes());

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleBiolEnvironmentTypes

  const ExampleClass::EnumType ExampleBiolEnvironmentTypes::s_Instance
  (
    GetExamples().AddEnum( ExampleBiolEnvironmentTypes())
  );

} // namespace bcl
