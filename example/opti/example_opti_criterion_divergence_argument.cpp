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
#include "opti/bcl_opti_criterion_divergence_argument.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_opti_divergence_argument.cpp
  //! @brief tests the implementation of CriterionDivergenceArgument
  //!
  //! @author fischea
  //! @date 12/13/2012
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleOptiCriterionDivergenceArgument :
    public ExampleInterface
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

  public:

    //! @brief Clone function
    //! @return pointer to a new ExampleOptiCriterionDivergenceArgument
    ExampleOptiCriterionDivergenceArgument *Clone() const
    {
      return new ExampleOptiCriterionDivergenceArgument( *this);
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
      storage::Pair< double, double>::s_Instance.IsDefined();
      opti::CriterionDivergenceArgument< double, double>::s_Instance.IsDefined();

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // create criterion with a maximum number of repeats of 3 and an interval width of 0.25
      util::ShPtr< opti::CriterionDivergenceArgument< double, double> > sp_criterion
      (
        new opti::CriterionDivergenceArgument< double, double>( 3, 0.25)
      );

      BCL_Example_Check
      (
        sp_criterion->GetMaxNumberRepeats() == 3 &&
        !util::IsDefined( sp_criterion->GetNumberRepeats()) &&
        sp_criterion->GetIntervalWidth() == 0.25,
        "members have not been set correctly. MaxNumberRepeats is " +
        util::Format()( sp_criterion->GetMaxNumberRepeats()) + " and should be 3." +
        " NumberRepeats is " + util::Format()( sp_criterion->GetNumberRepeats()) + " and should be 0." +
        " IntervalWidth is " + util::Format()( sp_criterion->GetIntervalWidth()) + " and should be 0.25."
      );

      // clone the criterion
      util::ShPtr< opti::CriterionDivergenceArgument< double, double> > sp_criterion_clone( sp_criterion->Clone());

      BCL_Example_Check
      (
        sp_criterion->GetMaxNumberRepeats() == sp_criterion_clone->GetMaxNumberRepeats() &&
        sp_criterion->GetNumberRepeats() == sp_criterion_clone->GetNumberRepeats() &&
        sp_criterion->GetIntervalWidth() == sp_criterion_clone->GetIntervalWidth(),
        "cloning the criterion was not successful"
      );

    /////////////////
    // data access //
    /////////////////

      BCL_ExampleCheck
      (
        sp_criterion->GetClassIdentifier(),
        ( GetStaticClassName< opti::CriterionDivergenceArgument< double, double> >())
      );

      BCL_Example_Check
      (
        sp_criterion->GetMaxNumberRepeats() == 3 &&
        !util::IsDefined( sp_criterion->GetNumberRepeats()) &&
        sp_criterion->GetIntervalWidth() == 0.25,
        "members have not been set correctly. MaxNumberRepeats is " +
        util::Format()( sp_criterion->GetMaxNumberRepeats()) + " and should be 3." +
        " NumberRepeats is " + util::Format()( sp_criterion->GetNumberRepeats()) + " and should be 0." +
        " IntervalWidth is " + util::Format()( sp_criterion->GetIntervalWidth()) + " and should be 0.25."
      );

    ////////////////
    // operations //
    ////////////////

      // create the tracker
      util::ShPtr< opti::Tracker< double, double> > sp_tracker
      (
        new opti::Tracker< double, double>( opti::e_SmallerIsBetter)
      );

      // create first tracking model
      const util::ShPtr< storage::Pair< double, double> > sp_model_first
      (
        new storage::Pair< double, double>( 1.0, 2.0)
      );

      // create second tracking model
      const util::ShPtr< storage::Pair< double, double> > sp_model_second
      (
        new storage::Pair< double, double>( 1.5, 2.0)
      );

      // create third tracking model
      const util::ShPtr< storage::Pair< double, double> > sp_model_third
      (
        new storage::Pair< double, double>( 1.1, 2.5)
      );

      // insert three models into the tracker and check if CriteriaMet is working
      sp_tracker->Track( sp_model_first);

      BCL_Example_Check
      (
        sp_criterion->CriteriaMet( *sp_tracker) == false,
        "CriteriaMet is true but should be false"
      );

      sp_tracker->Track( sp_model_first);

      BCL_Example_Check
      (
        sp_criterion->CriteriaMet( *sp_tracker) == false,
        "CriteriaMet is true but should be false"
      );

      sp_tracker->Track( sp_model_second);

      BCL_Example_Check
      (
        sp_criterion->CriteriaMet( *sp_tracker) == false,
        "CriteriaMet is true but should be false"
      );

      sp_tracker->Track( sp_model_second);

      BCL_Example_Check
      (
        sp_criterion->CriteriaMet( *sp_tracker) == false,
        "CriteriaMet is true but should be false"
      );

      sp_tracker->Track( sp_model_first);

      BCL_Example_Check
      (
        sp_criterion->CriteriaMet( *sp_tracker) == false,
        "CriteriaMet is true but should be false"
      );

      sp_tracker->Track( sp_model_second);

      BCL_Example_Check
      (
        sp_criterion->CriteriaMet( *sp_tracker) == false,
        "CriteriaMet is true but should be false"
      );

      sp_tracker->Track( sp_model_third);

      BCL_Example_Check
      (
        sp_criterion->CriteriaMet( *sp_tracker) == true,
        "CriteriaMet is false but should be true"
      );

    //////////////////////
    // input and output //
    //////////////////////

      // write criterion
      WriteBCLObject( *sp_criterion);

      // read criterion
      opti::CriterionDivergenceArgument< double, double> criterion_read;
      ReadBCLObject( criterion_read);

      // compare read in object to written object
      BCL_Example_Check
      (
        sp_criterion->GetMaxNumberRepeats() == criterion_read.GetMaxNumberRepeats() &&
        sp_criterion->GetNumberRepeats() == criterion_read.GetNumberRepeats() &&
        sp_criterion->GetIntervalWidth() == criterion_read.GetIntervalWidth(),
        "writing out and/or reading in the criterion was not successful"
      );

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; // class ExampleOptiCriterionDivergenceArgument

  const ExampleClass::EnumType ExampleOptiCriterionDivergenceArgument::s_Instance
  (
    GetExamples().AddEnum( ExampleOptiCriterionDivergenceArgument())
  );

} // namespace bcl
