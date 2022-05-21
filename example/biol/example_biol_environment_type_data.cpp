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
#include "biol/bcl_biol_environment_type_data.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_environment_types.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_biol_environment_type_data.cpp
  //!
  //! @author woetzen
  //! @date Nov 6, 2010 
  //! @remarks status complete
  //! @remarks reviewed by nobody on 
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleBiolEnvironmentTypeData :
    public ExampleInterface
  {
  public:

    ExampleBiolEnvironmentTypeData *Clone() const
    {
      return new ExampleBiolEnvironmentTypeData( *this);
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
      biol::EnvironmentTypeData env_type_data_default;

      BCL_Example_Check
      (
        env_type_data_default.GetTwoLetterCode() == ""
        && !util::IsDefined( env_type_data_default.GetReducedIndex())
        && !env_type_data_default.GetReducedType().IsDefined()
        && env_type_data_default.GetReducedTypeString() == ""
        && !util::IsDefined( env_type_data_default.GetDefaultThickness())
        && !env_type_data_default.IsGap(),
        "default constructor does not behave as expected"
      );

      // construct from information
      biol::EnvironmentTypeData env_type_data_a( "MC", "MC", 0, false, 10.0);

      // clone
      util::ShPtr< util::ObjectInterface> ptr( env_type_data_a.Clone());

    /////////////////
    // data access //
    /////////////////

      BCL_MessageStd
      (
        "static name of class is is " + GetStaticClassName< biol::EnvironmentTypeData>()
      );
      BCL_MessageStd( "class identifier is " + ptr->GetClassIdentifier());

      BCL_Example_Check
      (
        GetStaticClassName< biol::EnvironmentTypeData>() == ptr->GetClassIdentifier(),
        "incorrect class identifier"
      );

      BCL_MessageStd( "two letter code: " + env_type_data_a.GetTwoLetterCode());
      BCL_MessageStd( "reduced index: " + util::Format()( env_type_data_a.GetReducedIndex()));
      BCL_MessageStd( "reduced type: " + util::Format()( env_type_data_a.GetReducedType()));
      BCL_MessageStd( "reduced type string: " + util::Format()( env_type_data_a.GetReducedTypeString()));
      BCL_MessageStd( "default thickness: " + util::Format()( env_type_data_a.GetDefaultThickness()));
      BCL_MessageStd( "is gap: " + util::Format()( env_type_data_a.IsGap()));

      BCL_Example_Check
      (
        env_type_data_a.GetTwoLetterCode() == "MC"
        && env_type_data_a.GetReducedIndex() == 0
        && env_type_data_a.GetReducedType() == biol::GetEnvironmentTypes().e_MembraneCore
        && env_type_data_a.GetReducedTypeString() == biol::GetEnvironmentTypes().e_MembraneCore->GetTwoLetterCode()
        && env_type_data_a.GetDefaultThickness() == 10.0
        && !env_type_data_a.IsGap(),
        "constructor from information does not work correctly"
      );

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

      // write to file
      WriteBCLObject( env_type_data_a);
      biol::EnvironmentTypeData env_type_data_read;
      ReadBCLObject( env_type_data_read);

      BCL_Example_Check
      (
        env_type_data_a.GetTwoLetterCode() == env_type_data_read.GetTwoLetterCode()
        && env_type_data_a.GetReducedIndex() == env_type_data_read.GetReducedIndex()
        && env_type_data_a.GetReducedType() == env_type_data_read.GetReducedType()
        && env_type_data_a.GetReducedTypeString() == env_type_data_read.GetReducedTypeString()
        && env_type_data_a.GetDefaultThickness() == env_type_data_read.GetDefaultThickness()
        && env_type_data_a.IsGap() == env_type_data_read.IsGap(),
        "written and read env type data do not agree"
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleBiolEnvironmentTypeData

  const ExampleClass::EnumType ExampleBiolEnvironmentTypeData::s_Instance
  (
    GetExamples().AddEnum( ExampleBiolEnvironmentTypeData())
  );
  
} // namespace bcl
