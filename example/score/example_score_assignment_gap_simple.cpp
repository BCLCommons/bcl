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
#include "score/bcl_score_assignment_gap_simple.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_score_assignment_gap_simple.cpp
  //!
  //! @author alexanns
  //! @date
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on 
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleScoreAssignmentGapSimple :
    public ExampleInterface
  {
  public:

    ExampleScoreAssignmentGapSimple *Clone() const
    {
      return new ExampleScoreAssignmentGapSimple( *this);
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

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test default constructor
      BCL_MessageStd( "constructor from penalty");
      score::AssignmentGapSimple constructor_penalty( 0.0);
      BCL_Example_Check
      (
        constructor_penalty.GetPenalty() == 0, "constructor from penalty did not work properly"
      );

    /////////////////
    // data access //
    /////////////////

      // test SetData function and GetData functions
      BCL_MessageStd( "test SetData function and GetData functions");
      constructor_penalty.SetPenalty( 5.0);
      BCL_Example_Check
      (
        constructor_penalty.GetPenalty() == 5.0, "SetData function did not work properly"
      );

    ///////////////
    // operators //
    ///////////////

      // test operator()
      BCL_MessageStd( "test operator()");
      BCL_Example_Check
      (
        constructor_penalty( 10) == 50.0, "operator returned "
        + util::Format()( constructor_penalty( 10)) + " but should have return 50.0"
      );

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

      // write object
      WriteBCLObject( constructor_penalty);

      // read object
      score::AssignmentGapSimple gap_simple_read( util::GetUndefined< double>());
      ReadBCLObject( gap_simple_read);
      BCL_ExampleCheck( gap_simple_read.GetPenalty(), gap_simple_read.GetPenalty());

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleScoreAssignmentGapSimple

  const ExampleClass::EnumType ExampleScoreAssignmentGapSimple::s_Instance
  (
    GetExamples().AddEnum( ExampleScoreAssignmentGapSimple())
  );

} // namespace bcl
