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
#include "opti/bcl_opti_criterion_unimproved.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_opti_criterion_unimproved.cpp
  //! @brief tests the implementation of CriterionUnimproved
  //!
  //! @author fischea
  //! @date 1/4/2013
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleOptiCriterionUnimproved :
    public ExampleInterface
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

  public:

    //! @brief Clone function
    //! @return pointer to a new ExampleOptiCriterionUnimproved
    ExampleOptiCriterionUnimproved *Clone() const
    {
      return new ExampleOptiCriterionUnimproved( *this);
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

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // create criterion with maximum 3 unimproved steps
      opti::CriterionUnimproved< double, double> criterion( 4);

      BCL_ExampleIndirectCheck( criterion.GetMaxNumberUnimprovedSteps(), 4, "constructor");

      // create criterion with default constructor
      opti::CriterionUnimproved< double, double> criterion_default;

      BCL_ExampleIndirectCheck( criterion_default.GetMaxNumberUnimprovedSteps(), 0, "default constructor");

      // clone the criterion
      util::ShPtr< opti::CriterionUnimproved< double, double> > sp_criterion_clone( criterion.Clone());

      BCL_Example_Check
      (
        sp_criterion_clone->GetMaxNumberUnimprovedSteps() == 4,
        "cloning the criterion was not successful"
      );

    /////////////////
    // data access //
    /////////////////

      // check the retrieving of the class identifier
      BCL_ExampleCheck
      (
        criterion.GetClassIdentifier(),
        ( GetStaticClassName< opti::CriterionUnimproved< double, double> >())
      );

      // create criterion with other values
      opti::CriterionUnimproved< double, double> criterion_get( 5);

      // check the getter
      BCL_ExampleCheck( criterion_get.GetMaxNumberUnimprovedSteps(), 5);

    ////////////////
    // operations //
    ////////////////

      // create the tracker for the testing
      opti::Tracker< double, double> tracker( opti::e_SmallerIsBetter);

      for( size_t i( 0); i < 3; ++i)
      {
        util::ShPtr< storage::Pair< double, double> > sp_pair
        (
          new storage::Pair< double, double>( 1.0, 2.0)
        );
        tracker.Track( sp_pair);

        // check if criterion is met
        BCL_Example_Check
        (
          !criterion.CriteriaMet( tracker),
          "criterion is met but shouldn't"
        );
      }

      for( size_t i( 1); i < 6; ++i)
      {
        util::ShPtr< storage::Pair< double, double> > sp_pair
        (
          new storage::Pair< double, double>( 1.0, 1.0 + i * 0.1)
        );
        tracker.Track( sp_pair);

        // check if criterion is met
        BCL_Example_Check
        (
          !criterion.CriteriaMet( tracker),
          "criterion is met but shouldn't"
        );
      }

      for( size_t i( 0); i < 3; ++i)
      {
        util::ShPtr< storage::Pair< double, double> > sp_pair
        (
          new storage::Pair< double, double>( 1.0, 1.5 + i * 0.1)
        );
        tracker.Track( sp_pair);

        // check if criterion is met
        BCL_Example_Check
        (
          criterion.CriteriaMet( tracker),
          "criterion is not met but should be"
        );
      }

    //////////////////////
    // input and output //
    //////////////////////

      // write criterion
      WriteBCLObject( criterion);

      // read criterion
      opti::CriterionUnimproved< double, double> criterion_read;
      ReadBCLObject( criterion_read);

      // compare read in object to written object
      BCL_ExampleCheck( criterion_read.GetMaxNumberUnimprovedSteps(), criterion.GetMaxNumberUnimprovedSteps());

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; // class ExampleOptiCriterionUnimproved

  const ExampleClass::EnumType ExampleOptiCriterionUnimproved::s_Instance
  (
    GetExamples().AddEnum( ExampleOptiCriterionUnimproved())
  );

} // namespace bcl
