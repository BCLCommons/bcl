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
#include "opti/bcl_opti_criterion_n_step.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_opti_criterion_n_step.cpp
  //! @brief tests the implementation of CriterionNStep
  //!
  //! @author fischea
  //! @date Dec 12, 2012
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleOptiCriterionNStep :
    public ExampleInterface
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

  public:

    //! @brief Clone function
    //! @return pointer to a new ExampleOptiCriterionNStep
    ExampleOptiCriterionNStep *Clone() const
    {
      return new ExampleOptiCriterionNStep( *this);
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

    /////////////////
    // preparation //
    /////////////////

      // create the tracker
      opti::Tracker< int, int> tracker;

      // create the tracking objectove for the testing
      const util::ShPtr< storage::Pair< int, int> > sp_pair( new storage::Pair< int, int>( 1, 1));

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // create the criterion with a repeat interval of two
      opti::CriterionNStep< int, int> criterion( 2);

    /////////////////
    // data access //
    /////////////////

      // check if the correct repeat interval is returned
      BCL_ExampleCheck( criterion.GetRepeatInterval(), 2);

      // test getter for class name identifier
      BCL_ExampleCheck
      (
        criterion.GetClassIdentifier(),
        ( GetStaticClassName< opti::CriterionNStep< int, int> >())
      );

    ////////////////
    // operations //
    ////////////////

      tracker.Track( sp_pair);
      BCL_ExampleCheck( criterion.CriteriaMet( tracker), false);

      tracker.Track( sp_pair);
      BCL_ExampleCheck( criterion.CriteriaMet( tracker), true);

      tracker.Track( sp_pair);
      BCL_ExampleCheck( criterion.CriteriaMet( tracker), false);

      tracker.Track( sp_pair);
      BCL_ExampleCheck( criterion.CriteriaMet( tracker), true);

    //////////////////////
    // input and output //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; // class ExampleOptiCriterionNStep

  const ExampleClass::EnumType ExampleOptiCriterionNStep::s_Instance
  (
    GetExamples().AddEnum( ExampleOptiCriterionNStep())
  );

} // namespace bcl
