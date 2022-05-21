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
#include "opti/bcl_opti_criterion_phase.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_opti_criterion_phase.cpp
  //! @brief tests the implementation of CriterionPhase
  //!
  //! @author fischea
  //! @date Dec 14, 2013
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleOptiCriterionPhase :
    public ExampleInterface
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

  public:

    //! @brief Clone function
    //! @return pointer to a new ExampleOptiCriterionPhase
    ExampleOptiCriterionPhase *Clone() const
    {
      return new ExampleOptiCriterionPhase( *this);
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

    //////////////////
    // preparations //
    //////////////////

      // create sets that contains the phases for the criterion
      storage::Set< opti::PhaseEnum> phases_end;
      phases_end.InsertElement( opti::e_End);
      storage::Set< opti::PhaseEnum> phases_always;
      phases_always.InsertElement( opti::e_Always);

      // create the tracker
      opti::Tracker< int, int> tracker;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // construct the criteria
      opti::CriterionPhase< int, int> criterion_end( phases_end);
      opti::CriterionPhase< int, int> criterion_always( phases_always);

    /////////////////
    // data access //
    /////////////////

      // test getter for class name identifier
      BCL_ExampleCheck
      (
        criterion_end.GetClassIdentifier(), ( GetStaticClassName< opti::CriterionPhase< int, int> >())
      );

    ////////////////
    // operations //
    ////////////////

      tracker.SetPhase( opti::e_Start);
      BCL_ExampleCheck( criterion_end.CriteriaMet( tracker), false);
      BCL_ExampleCheck( criterion_always.CriteriaMet( tracker), true);

      tracker.SetPhase( opti::e_Iteration);
      BCL_ExampleCheck( criterion_end.CriteriaMet( tracker), false);
      BCL_ExampleCheck( criterion_always.CriteriaMet( tracker), true);

      tracker.SetPhase( opti::e_End);
      BCL_ExampleCheck( criterion_end.CriteriaMet( tracker), true);
      BCL_ExampleCheck( criterion_always.CriteriaMet( tracker), true);

    //////////////////////
    // input and output //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; // class ExampleOptiCriterionPhase

  const ExampleClass::EnumType ExampleOptiCriterionPhase::s_Instance
  (
    GetExamples().AddEnum( ExampleOptiCriterionPhase())
  );

} // namespace bcl
