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
#include "opti/bcl_opti_criterion_improvement_ratio.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_opti_criterion_improvement_ratio.cpp
  //! @brief tests the implementation of CriterionImprovementRatio
  //!
  //! @author mendenjl
  //! @date Sep 03, 2013
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleOptiCriterionImprovementRatio :
    public ExampleInterface
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

  public:

    //! @brief Clone function
    //! @return pointer to a new ExampleOptiCriterionImprovementRatio
    ExampleOptiCriterionImprovementRatio *Clone() const
    {
      return new ExampleOptiCriterionImprovementRatio( *this);
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
//      // construct empty nodes list
//      const util::SiPtrList< const math::MutateInterface< double > > nodes;
//
//      // construct starting argument and step
//      util::ShPtr< storage::Pair< util::ShPtr< double>, double>> start_argument( new storage::Pair< util::ShPtr< double>, double>( 100, 0));
//      mc::Step< double, double> start_step( start_argument, 0, opti::e_Improved, nodes);
//
//      // construct arguments
//      util::ShPtr< storage::Pair< util::ShPtr< double>, double> > arg_1( new storage::Pair< util::ShPtr< double>, double>( 101, -10));
//      util::ShPtr< storage::Pair< util::ShPtr< double>, double> > arg_2( new storage::Pair< util::ShPtr< double>, double>( 102, -11));
//      util::ShPtr< storage::Pair< util::ShPtr< double>, double> > arg_3( new storage::Pair< util::ShPtr< double>, double>( 103, -15));
//      util::ShPtr< storage::Pair< util::ShPtr< double>, double> > arg_4( new storage::Pair< util::ShPtr< double>, double>( 104,  -8));
//      util::ShPtr< storage::Pair< util::ShPtr< double>, double> > arg_5( new storage::Pair< util::ShPtr< double>, double>( 105, -17));
//      util::ShPtr< storage::Pair< util::ShPtr< double>, double> > arg_6( new storage::Pair< util::ShPtr< double>, double>( 106, -13));
//      util::ShPtrVector< storage::Pair< util::ShPtr< double>, double>> arguments;
//      arguments.PushBack( arg_1);
//      arguments.PushBack( arg_2);
//      arguments.PushBack( arg_3);
//      arguments.PushBack( arg_4);
//      arguments.PushBack( arg_5);
//      arguments.PushBack( arg_6);
//
//      // construct steps
//      storage::Vector< mc::Step< double, double> > steps;
//      steps.PushBack( mc::Step< double, double>( arg_1, 1, opti::e_Improved, nodes));
//      steps.PushBack( mc::Step< double, double>( arg_2, 2, opti::e_Improved, nodes));
//      steps.PushBack( mc::Step< double, double>( arg_3, 3, opti::e_Improved, nodes));
//      steps.PushBack( mc::Step< double, double>( arg_4, 4, opti::e_Rejected, nodes));
//      steps.PushBack( mc::Step< double, double>( arg_5, 5, opti::e_Improved, nodes));
//      steps.PushBack( mc::Step< double, double>( arg_6, 6, opti::e_Rejected, nodes));
//
//      // construct expect results
//      storage::Vector< size_t> expected_criteria_met;
//      expected_criteria_met.PushBack( 0);
//      expected_criteria_met.PushBack( 0);
//      expected_criteria_met.PushBack( 0);
//      expected_criteria_met.PushBack( 0);
//      expected_criteria_met.PushBack( 0);
//      expected_criteria_met.PushBack( 1);
//
//      // construct expected last satisfactory result
//      storage::Vector< double> expected_last_satisfactory_results;
//      expected_last_satisfactory_results.PushBack( -10);
//      expected_last_satisfactory_results.PushBack( -10);
//      expected_last_satisfactory_results.PushBack( -15);
//      expected_last_satisfactory_results.PushBack( -15);
//      expected_last_satisfactory_results.PushBack( -15);
//      expected_last_satisfactory_results.PushBack( -15);
//
//      // construct expected last satisfactory step number
//      storage::Vector< size_t> expected_last_satisfactory_step_numbers;
//      expected_last_satisfactory_step_numbers.PushBack( 1);
//      expected_last_satisfactory_step_numbers.PushBack( 1);
//      expected_last_satisfactory_step_numbers.PushBack( 3);
//      expected_last_satisfactory_step_numbers.PushBack( 3);
//      expected_last_satisfactory_step_numbers.PushBack( 3);
//      expected_last_satisfactory_step_numbers.PushBack( 3);
//
//      // construct a tracker
//      mc::Tracker< double, double> tracker( 1);
//
//      // insert the starting argument into the tracker
//      tracker.Insert( start_argument);
//      tracker.InsertMCStep( start_step);
//
//    //////////////////////////////////
//    // construction and destruction //
//    //////////////////////////////////
//
//      // ratio
//      const double ratio( 0.2);
//      const size_t nr_steps( 3);
//
//      // default constructor
//      BCL_MessageStd( "test default constructor");
//      mc::TerminateImprovementRatio< double, double> terminate_def;
//
//      // constructor from boolean and strand counts
//      BCL_MessageStd( "test constructor");
//      mc::TerminateImprovementRatio< double, double> terminate( ratio, nr_steps, tracker);
//      BCL_Example_Check
//      (
//        terminate.GetRatio() == ratio &&
//        terminate.GetNumberSteps() == nr_steps &&
//        terminate.GetTracker().GetPointer() == &tracker &&
//        terminate.GetLastSatisfactoryResult() == 0 &&
//        terminate.GetLastSatisfactoryStepNumber() == 0,
//        "The default constructor failed " + util::Format()( terminate)
//      );
//
//      BCL_MessageStd( "test clone constructor");
//      util::ShPtr< mc::TerminateImprovementRatio< double, double> > sp_terminate( terminate.Clone());
//      BCL_Example_Check
//      (
//        sp_terminate->GetRatio() == terminate.GetRatio() &&
//        sp_terminate->GetNumberSteps() == terminate.GetNumberSteps() &&
//        sp_terminate->GetTracker().GetPointer() == terminate.GetTracker().GetPointer() &&
//        sp_terminate->GetLastSatisfactoryResult() == terminate.GetLastSatisfactoryResult() &&
//        sp_terminate->GetLastSatisfactoryStepNumber() == terminate.GetLastSatisfactoryStepNumber(),
//        "The cloned mutate is different!\n" + util::Format()( *sp_terminate) + "\nvs\n" + util::Format()( terminate)
//      );
//
//    /////////////////
//    // data access //
//    /////////////////
//
//      // test the GetClassIdentifier
//      BCL_ExampleCheck( terminate.GetClassIdentifier(), "bcl::mc::TerminateImprovementRatio<double,double>");
//
//      // test GetRatio()
//      BCL_ExampleCheck( terminate.GetRatio(), ratio);
//
//      // test GetTracker()
//      BCL_MessageStd( "test GetTracker()");
//      BCL_Example_Check
//      (
//        terminate.GetTracker().GetPointer() == &tracker,
//        "The GetTracker() should return a pointer of adress " + util::Format()( &tracker) +
//        " not " + util::Format()( terminate.GetTracker().GetPointer())
//      );
//
//      // test GetLastSatisfactoryResult()
//      BCL_ExampleCheck( terminate.GetLastSatisfactoryResult(), 0);
//
//      // test GetLastSatisfactoryStepNumber()
//      BCL_ExampleCheck( terminate.GetLastSatisfactoryStepNumber(), 0);
//
//    ///////////////
//    // operators //
//    ///////////////
//
//      // initialize vector to hold the results back
//      storage::Vector< size_t> criteria_met;
//      storage::Vector< double> last_satisfactory_results;
//      storage::Vector< size_t> last_satisfactory_step_numbers;
//
//      BCL_MessageStd( "test CriteriaMet() in an iterative fashion");
//      for( size_t i( 0); i < steps.GetSize(); ++i)
//      {
//        // insert the step into tracker
//        BCL_MessageStd( "inserting step #" + util::Format()( i + i) + " into the tracker");
//        tracker.InsertMCStep( steps( i));
//
//        // if this was an improved step insert into opti::Tracker
//        if( steps( i).GetStatus() == opti::e_Improved)
//        {
//          tracker.Insert( arguments( i));
//        }
//
//        // store variables
//        const bool this_status( terminate.CriteriaMet());
//        const double this_result( terminate.GetLastSatisfactoryResult());
//        const size_t this_step_number( terminate.GetLastSatisfactoryStepNumber());
//
//        // output values
//        BCL_MessageStd
//        (
//          "criteria_met: " + util::Format()( this_status) +
//          " last result: " + util::Format()( this_result) +
//          " last_step_number: " + util::Format()( this_step_number)
//        );
//
//        // insert the results
//        criteria_met.PushBack( this_status);
//        last_satisfactory_results.PushBack( this_result);
//        last_satisfactory_step_numbers.PushBack( this_step_number);
//      }
//
//      // now check the expected values with the result
//      BCL_Example_Check
//      (
//        criteria_met == expected_criteria_met,
//        "The CriteriaMet() results are different from expected\n " +
//        util::Format()( criteria_met) + "\nvs\n" + util::Format()( expected_criteria_met)
//      );
//
//      BCL_Example_Check
//      (
//        last_satisfactory_results == expected_last_satisfactory_results,
//        "The last satisfactory results are different from expected\n " +
//        util::Format()( last_satisfactory_results) + "\nvs\n" + util::Format()( expected_last_satisfactory_results)
//      );
//
//      BCL_Example_Check
//      (
//        last_satisfactory_step_numbers == expected_last_satisfactory_step_numbers,
//        "The last satisfactory step numbers are different from expected\n " +
//        util::Format()( last_satisfactory_step_numbers) + "\nvs\n" + util::Format()( expected_last_satisfactory_step_numbers)
//      );
//
//    ////////////////
//    // operations //
//    ////////////////
//
//    //////////////////////
//    // input and output //
//    //////////////////////
//
//      // test read write
//      BCL_MessageStd( "Testing read write")
//      WriteBCLObject( terminate);
//      mc::TerminateImprovementRatio< double, double> terminate_read;
//      ReadBCLObject( terminate_read);
//      BCL_Example_Check
//      (
//        terminate.GetRatio() == terminate_read.GetRatio() &&
//        terminate.GetNumberSteps() == terminate_read.GetNumberSteps(),
//        "The read terminate is different\n" + util::Format()( terminate_read) + "\nvs\n" + util::Format()( terminate)
//      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

   }; // class ExampleOptiCriterionImprovementRatio

   const ExampleClass::EnumType ExampleOptiCriterionImprovementRatio::s_Instance
   (
     GetExamples().AddEnum( ExampleOptiCriterionImprovementRatio())
   );

} // namespace bcl
