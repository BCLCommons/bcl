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
#include "opti/bcl_opti_criterion_result_threshold.h"

// includes from bcl - sorted alphabetically
#include "opti/bcl_opti_tracker.h"

// external includes - sorted alphabetically

namespace bcl
{
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_opti_criterion_result_threshold.cpp
  //! @brief tests the implementation of the CriterionResultThreshold
  //!
  //! @author mendenjl
  //! @date Aug 30, 2013
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleOptiCriterionResultThreshold :
    public ExampleInterface
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

  public:

    //! @brief Clone function
    //! @return pointer to a new ExampleOptiCriterionResultThreshold
    ExampleOptiCriterionResultThreshold *Clone() const
    {
      return new ExampleOptiCriterionResultThreshold( *this);
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
      // create a tracker
      opti::Tracker< double, double> tracker( opti::e_SmallerIsBetter);

      // track models for testing
      tracker.Track
      (
        util::ShPtr< storage::Pair< double, double> >
        (
          new storage::Pair< double, double>( 3.0, 3.0)
        )
      );

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      opti::CriterionResultThreshold< double, double> criterion_default;

      // construct from threshold
      opti::CriterionResultThreshold< double, double> criterion( 0.5);

      // clone function
      util::ShPtr< opti::CriterionResultThreshold< double, double> > criterion_clone
      (
        criterion.Clone()
      );

    /////////////////
    // data access //
    /////////////////

      BCL_ExampleCheck
      (
        criterion.GetClassIdentifier(),
        ( GetStaticClassName< opti::CriterionResultThreshold< double, double> >())
      );

    ////////////////
    // operations //
    ////////////////

      BCL_Example_Check
      (
        criterion.CriteriaMet( tracker) == false,
        "incorrect termination: should be false, but is true"
      );

      // insert new argument result pair into the tracker
      tracker.Track
      (
        util::ShPtr< storage::Pair< double, double> >
        (
          new storage::Pair< double, double>( 3.0, 0.4)
        )
      );

      BCL_Example_Check
      (
        criterion.CriteriaMet( tracker) == true,
        "incorrect termination: should be true, but is false"
      );

    //////////////////////
    // input and output //
    //////////////////////

      // write criterion
      WriteBCLObject( criterion);

      // read in criterion
      opti::CriterionResultThreshold< double, double> criterion_read;
      ReadBCLObject( criterion_read);

      // check class identifier
      BCL_Example_Check
      (
        criterion_read.CriteriaMet( tracker) == true,
        "writing out or reading in of criterion failed"
      );

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; // class ExampleOptiCriterionResultThreshold

  const ExampleClass::EnumType ExampleOptiCriterionResultThreshold::s_Instance
  (
    GetExamples().AddEnum( ExampleOptiCriterionResultThreshold())
  );

} // namespace bcl
