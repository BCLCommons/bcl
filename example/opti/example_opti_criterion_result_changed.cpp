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
#include "opti/bcl_opti_criterion_result_changed.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_opti_criterion_result_changed.cpp
  //! @brief tests the implementation of CriterionResultChanged
  //!
  //! @author mendenjl
  //! @date Jul 18, 2017
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleOptiCriterionResultChanged :
    public ExampleInterface
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

  public:

    //! @brief Clone function
    //! @return pointer to a new ExampleOptiCriterionResultChanged
    ExampleOptiCriterionResultChanged *Clone() const
    {
      return new ExampleOptiCriterionResultChanged( *this);
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
      opti::CriterionResultChanged< double, double> criterion;

      // clone the criterion
      util::ShPtr< opti::CriterionResultChanged< double, double> > sp_criterion_clone( criterion.Clone());

    /////////////////
    // data access //
    /////////////////

      // check the retrieving of the class identifier
      BCL_ExampleCheck
      (
        criterion.GetClassIdentifier(),
        ( GetStaticClassName< opti::CriterionResultChanged< double, double> >())
      );

    ////////////////
    // operations //
    ////////////////

      // create the tracker for the testing
      opti::Tracker< double, double> tracker( opti::e_SmallerIsBetter);
      util::ShPtr< storage::Pair< double, double> > sp_pair
      (
        new storage::Pair< double, double>( 1.0, 2.0)
      );
      tracker.Track( sp_pair);
      BCL_ExampleIndirectCheck
      (
        criterion.CriteriaMet( tracker), true,
        "first insert should trigger that result changed"
      );
      for( size_t i( 0); i < 3; ++i)
      {
        tracker.Track( sp_pair);

        // check if criterion is met
        BCL_ExampleIndirectCheck
        (
          criterion.CriteriaMet( tracker),
          false,
          "adding back same result should not trigger (result unchanged)"
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
        BCL_ExampleIndirectCheck
        (
          criterion.CriteriaMet( tracker),
          true,
          "adding back differing result should trigger"
        );
      }

    //////////////////////
    // input and output //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; // class ExampleOptiCriterionResultChanged

  const ExampleClass::EnumType ExampleOptiCriterionResultChanged::s_Instance
  (
    GetExamples().AddEnum( ExampleOptiCriterionResultChanged())
  );

} // namespace bcl
