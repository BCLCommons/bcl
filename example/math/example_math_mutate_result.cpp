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
#include "math/bcl_math_mutate_result.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_math_mutate_result.cpp
  //!
  //! @author karakam
  //! @date
  //! @remarks status empty
  //! @remarks reviewed by nobody on 
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleMathMutateResult :
    public ExampleInterface
  {
  public:

    ExampleMathMutateResult *Clone() const
    {
      return new ExampleMathMutateResult( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! returns class name
    //! the class name as const ref std::string
    const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    int Run() const
    {

      // initialize an amino acid
      biol::AA ala_1( util::ShPtr< biol::AAData>( new biol::AAData( biol::GetAATypes().ALA, 1, 1)));
      biol::AA val_2( util::ShPtr< biol::AAData>( new biol::AAData( biol::GetAATypes().VAL, 2, 2)));

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      BCL_MessageStd( "Testing default constructor");
      math::MutateResult< biol::AA> default_result;

    /////////////////
    // data access //
    /////////////////

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
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleMathMutateResult

  const ExampleClass::EnumType ExampleMathMutateResult::s_Instance
  (
    GetExamples().AddEnum( ExampleMathMutateResult())
  );

} // namespace bcl

