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
#include "opti/bcl_opti_evolution_population_member.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_opti_evolution_population_member.cpp
  //! @brief this example tests the implementation of EvolutionPopulationMember
  //!
  //! @author geanesar
  //! @date Nov 1, 2014
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleOptiEvolutionPopulationMember :
     public ExampleInterface
   {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

  public:

    //! @brief Clone function
    //! @return pointer to a new ExampleOptiEvolutionPopulationMember
    ExampleOptiEvolutionPopulationMember *Clone() const
    {
      return new ExampleOptiEvolutionPopulationMember( *this);
    }

  //////////
  // data //
  //////////

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
      size_t member = 5, fitness = 7;
      util::ObjectDataLabel history( "test");

      // Construction
      opti::EvolutionPopulationMember< size_t, size_t> pop_member_def;
      opti::EvolutionPopulationMember< size_t, size_t> pop_member_member( member, fitness);
      opti::EvolutionPopulationMember< size_t, size_t> pop_member_full( member, fitness, history);

      // Set up everything
      pop_member_def.SetMember( member);
      pop_member_def.SetFitness( fitness);
      pop_member_def.SetHistory( history);
      pop_member_member.SetHistory( history);

      // Check that member info was set correctly
      BCL_ExampleCheck( pop_member_member.GetMember(), member);
      BCL_ExampleCheck( pop_member_full.GetMember(), member);
      BCL_ExampleCheck( pop_member_def.GetMember(), member);

      // Check that fitnesses were set correctly
      BCL_ExampleCheck( pop_member_def.GetFitness(), fitness);
      BCL_ExampleCheck( pop_member_member.GetFitness(), fitness);
      BCL_ExampleCheck( pop_member_full.GetFitness(), fitness);

      // Check that histories were set correctly
      BCL_ExampleCheck( pop_member_def.GetHistory(), history);
      BCL_ExampleCheck( pop_member_member.GetHistory(), history);
      BCL_ExampleCheck( pop_member_full.GetHistory(), history);

      // Check cloning
      util::ShPtr< opti::EvolutionPopulationMember< size_t, size_t> > new_member( pop_member_def.Clone());

      BCL_ExampleCheck
      (
        pop_member_def.GetMember() == new_member->GetMember()
        && pop_member_def.GetFitness() == new_member->GetFitness()
        && pop_member_def.GetHistory() == new_member->GetHistory(),
        true
      );

      pop_member_full.SetFitness( fitness - 1);

      // Check that setting the fitness works
      BCL_ExampleCheck( pop_member_full.GetFitness(), fitness - 1);

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

   }; // class ExampleOptiEvolutionPopulationMember

   const ExampleClass::EnumType ExampleOptiEvolutionPopulationMember::s_Instance
   (
     GetExamples().AddEnum( ExampleOptiEvolutionPopulationMember())
   );

} // namespace bcl
