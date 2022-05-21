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
#include "opti/bcl_opti_criterion_all.h"

// includes from bcl - sorted alphabetically
#include "opti/bcl_opti_criterion_convergence_result.h"
#include "opti/bcl_opti_criterion_number_iterations.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_opti_criterion_all.cpp
  //! @brief this example tests the implementation of the CriterionAll
  //!
  //! @author mendenjl
  //! @date Jul 18, 2017
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleOptiCriterionAll :
    public ExampleInterface
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

  public:

    //! @brief Clone function
    //! @return pointer to a new ExampleOptiCriterionAll
    ExampleOptiCriterionAll *Clone() const
    {
      return new ExampleOptiCriterionAll( *this);
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
      opti::CriterionAll< double, double>::s_Instance.IsDefined();
      opti::CriterionAll< double, double>::s_Instance.IsDefined();
      opti::CriterionNumberIterations< double, double>::s_Instance.IsDefined();

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // maximum number of iterations
      const size_t max_num_iterations( 5);

      // create criterion with a maximum number of iterations
      util::Implementation< opti::CriterionInterface< double, double> > sp_criterion_num_iterations
      (
        new opti::CriterionNumberIterations< double, double>( max_num_iterations)
      );

      // maximum number of repeats
      const size_t max_num_repeats( 3);

      // tolerance for the repeats
      const double tolerance( 0.25);

      // create criterion that terminate upon convergence of the result
      util::Implementation< opti::CriterionInterface< double, double> > sp_criterion_convergence_result
      (
        new opti::CriterionConvergenceResult< double, double>( max_num_repeats, tolerance)
      );

      // create combined criterion and add criteria to it
      storage::List< util::Implementation< opti::CriterionInterface< double, double> > > criteria_list;
      criteria_list.PushBack( sp_criterion_num_iterations);
      criteria_list.PushBack( sp_criterion_convergence_result);
      util::ShPtr< opti::CriterionAll< double, double> > sp_criterion_combine
      (
        new opti::CriterionAll< double, double>( criteria_list)
      );

    /////////////////
    // data access //
    /////////////////

      BCL_ExampleCheck
      (
        sp_criterion_combine->GetClassIdentifier(),
        ( GetStaticClassName< opti::CriterionAll< double, double> >())
      );

    ////////////////
    // operations //
    ////////////////

      // create the tracker
      util::ShPtr< opti::Tracker< double, double> > sp_tracker
      (
        new opti::Tracker< double, double>( opti::e_SmallerIsBetter)
      );

      // create tracking model
      const util::ShPtr< storage::Pair< double, double> > sp_model_first
      (
        new storage::Pair< double, double>( 1.0, 2.0)
      );

      // track two models and check if CriteriaMet is working
      for( size_t model_number( 0); model_number < 2; ++model_number)
      {
        sp_tracker->Track( sp_model_first);

        BCL_Example_Check
        (
          sp_criterion_combine->CriteriaMet( *sp_tracker) == false,
          "CriteriaMet is true but should be false"
        );
      }

      // create tracking model
      const util::ShPtr< storage::Pair< double, double> > sp_model_second
      (
        new storage::Pair< double, double>( 1.0, 3.0)
      );

      // track two models and check if CriteriaMet is working
      for( size_t model_number( 0); model_number < 5; ++model_number)
      {
        sp_tracker->Track( sp_model_second);

        BCL_Example_Check
        (
          sp_criterion_combine->CriteriaMet( *sp_tracker) == false,
          "CriteriaMet is true but should be false"
        );
      }

      // track again to hit the iteration limit and the convergence limit
      sp_tracker->Track( sp_model_second);
      BCL_Example_Check
      (
        sp_criterion_combine->CriteriaMet( *sp_tracker) == true,
        "CriteriaMet is false but should be true"
      );

    //////////////////////
    // input and output //
    //////////////////////

      // write criterion
      WriteBCLObject( *sp_criterion_combine);

      // read criterion
      opti::CriterionAll< double, double> criterion_read;
      ReadBCLObject( criterion_read);

      BCL_Example_Check
      (
        sp_criterion_combine->CriteriaMet( *sp_tracker) == criterion_read.CriteriaMet( *sp_tracker),
        "read in object does not match the written out object"
      );

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; // class ExampleOptiCriterionAll

  const ExampleClass::EnumType ExampleOptiCriterionAll::s_Instance
  (
    GetExamples().AddEnum( ExampleOptiCriterionAll())
  );

} // namespace bcl
