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
#include "opti/bcl_opti_criterion_number_iterations.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_opti_criterion_number_iterations.cpp
  //! @brief tests the implementation of CriterionNumberIterations
  //!
  //! @author fischea
  //! @date Dec 12, 2012
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleOptiCriterionNumberIterations :
    public ExampleInterface
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

  public:

    //! @brief Clone function
    //! @return pointer to a new ExampleOptiCriterionNumberIterations
    ExampleOptiCriterionNumberIterations *Clone() const
    {
      return new ExampleOptiCriterionNumberIterations( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns the class name
    //! @return the class name as const ref std::string
    const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief run routine
    //! this is performing the execution of the example
    int Run() const
    {
      opti::CriterionNumberIterations< int, int>::s_Instance.IsDefined();

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // create criterion with a maximum number of iterations of 10
      util::ShPtr< opti::CriterionNumberIterations< int, int> > sp_criterion
      (
        new opti::CriterionNumberIterations< int, int>( 10)
      );

      BCL_ExampleCheck( sp_criterion->GetMaxNumberIterations(), 10);

      // clone the criterion
      util::ShPtr< opti::CriterionNumberIterations< int, int> > sp_criterion_clone( sp_criterion->Clone());

      BCL_ExampleCheck( sp_criterion_clone->GetMaxNumberIterations(), 10);

    /////////////////
    // data access //
    /////////////////

      // set maximum number of iterations to 5
      sp_criterion->SetMaxNumberIterations( 5);

      BCL_ExampleCheck( sp_criterion->GetMaxNumberIterations(), 5);

      // test getter for class name identifier
      BCL_ExampleCheck
      (
        sp_criterion->GetClassIdentifier(),
        ( GetStaticClassName< opti::CriterionNumberIterations< int, int> >())
      );

    ////////////////
    // operations //
    ////////////////

      // create the tracker
      util::ShPtr< opti::Tracker< int, int> > sp_tracker
      (
        new opti::Tracker< int, int>( opti::e_SmallerIsBetter)
      );

      // create tracking model
      const util::ShPtr< storage::Pair< int, int> > sp_model( new storage::Pair< int, int>( 1, 2));

      // insert four models into the tracker and check if CriteriaMet is working
      for( size_t model_number( 0); model_number < 4; ++model_number)
      {
        sp_tracker->Track( sp_model);

        BCL_ExampleCheck( sp_criterion->CriteriaMet( *sp_tracker), false);
      }

      // insert more models into the tracker and check if CriteriaMet is working
      for( size_t model_number( 0); model_number < 4; ++model_number)
      {
        sp_tracker->Track( sp_model);

        BCL_ExampleCheck( sp_criterion->CriteriaMet( *sp_tracker), true);
      }

    //////////////////////
    // input and output //
    //////////////////////

      // write criterion
      WriteBCLObject( *sp_criterion);

      // read criterion
      opti::CriterionNumberIterations< int, int> criterion_read;
      ReadBCLObject( criterion_read);

      // compare read in object to written object
      BCL_ExampleCheck( criterion_read.GetMaxNumberIterations(), sp_criterion->GetMaxNumberIterations());

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; // class ExampleOptiCriterionNumberIterations

  const ExampleClass::EnumType ExampleOptiCriterionNumberIterations::s_Instance
  (
    GetExamples().AddEnum( ExampleOptiCriterionNumberIterations())
  );

} // namespace bcl
