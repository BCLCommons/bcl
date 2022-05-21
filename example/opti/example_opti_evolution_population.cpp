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
#include "opti/bcl_opti_evolution_population.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_opti_evolution_population.cpp
  //! @brief this example tests the evolution population class for functionality
  //!
  //! @author geanesar
  //! @date Nov 1 2014
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleOptiEvolutionPopulation :
     public ExampleInterface
   {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

  public:

    //! @brief Clone function
    //! @return pointer to a new ExampleOptiEvolutionPopulation
    ExampleOptiEvolutionPopulation *Clone() const
    {
      return new ExampleOptiEvolutionPopulation( *this);
    }

  //////////
  // data //
  //////////

    //! @brief Class that tests for uniqueness
    class SizeTPopulationUnique :
      public opti::EvolutionMemberUniqueInterface< size_t, size_t>
    {
    public:

        //! @brief default constructor
        SizeTPopulationUnique()
        {
        }

        //! @brief Clone function
        //! @return pointer to a new ExampleOptiApproximatorEvolution
        SizeTPopulationUnique *Clone() const
        {
          return new SizeTPopulationUnique( *this);
        }

        //! @brief returns the class name
        //! @return the class name as const ref std::string
        const std::string &GetClassIdentifier() const
        {
          return GetStaticClassName( *this);
        }

        //! @brief function for testing if a member is contained in the population
        //! @param MEMBER the member to check
        //! @param POPULATIon the population to search
        //! @return true if the member is contained in the population, false otherwise
        bool Contains
        (
          const opti::EvolutionPopulationMember< size_t, size_t> &MEMBER,
          const opti::EvolutionPopulation< size_t, size_t> &POPULATION
        ) const
        {
          iterate::Generic< const opti::EvolutionPopulationMember< size_t, size_t> > itr_member( POPULATION.GetMembersIterator());
          size_t member_value( MEMBER.GetMember());
          for( ; itr_member.NotAtEnd(); ++itr_member)
          {
            if( itr_member->GetMember() == member_value)
            {
              return true;
            }
          }
          return false;
        }

        //! @brief return parameters for member data that are set up from the labels
        //! @return parameters for member data that are set up from the labels
        io::Serializer GetSerializer() const
        {
          return io::Serializer();
        }

    };

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

      // Construction and destruction
      opti::EvolutionPopulation< size_t, size_t> pop;
      
      opti::EvolutionPopulation< size_t, size_t> pop_assign;
      pop_assign = pop;

      util::ShPtr< opti::EvolutionPopulation< size_t, size_t> > pop_copy
      ( 
        new opti::EvolutionPopulation< size_t, size_t>( pop)
      );

      pop_copy = util::ShPtr< opti::EvolutionPopulation< size_t, size_t> >();

      opti::EvolutionPopulationMember< size_t, size_t> test( 0, 0);
      util::ShPtrVector< opti::EvolutionPopulationMember< size_t, size_t> > init_members;
      init_members.PushBack
      ( 
        util::ShPtr< opti::EvolutionPopulationMember< size_t, size_t> >
        ( 
          new opti::EvolutionPopulationMember< size_t, size_t>( test)
        )
      ); 
      init_members.PushBack
      ( 
        util::ShPtr< opti::EvolutionPopulationMember< size_t, size_t> >
        ( 
          new opti::EvolutionPopulationMember< size_t, size_t>( test)
        )
      ); 

      opti::EvolutionPopulation< size_t, size_t> pop_init( 2, init_members);
      BCL_ExampleCheck( pop_init.GetCurrentSize(), 2);

      // Input/output
      BCL_ExampleCheck( pop.GetClassIdentifier(), "bcl::opti::EvolutionPopulation<size_t,size_t>");

      // Size setting/querying
      BCL_ExampleCheck( pop.GetNormalSize(), 0);

      pop.SetNormalSize( 1);
      BCL_ExampleCheck( pop.GetNormalSize(), 1);

      // Resetting should not change the normal size, but should remove members
      pop.Reset();
      BCL_ExampleCheck( pop.GetCurrentSize(), 0);
      BCL_ExampleCheck( pop.GetNormalSize(), 1);

      // Member addition
      pop.AddMember( size_t( 0), size_t( 1));
      BCL_ExampleCheck( pop.GetCurrentSize(), 1);

      pop.AddMember( size_t( 1), size_t( 1));
      BCL_ExampleCheck( pop.GetHighestFitness(), 1);
      BCL_ExampleCheck( pop.GetLowestFitness(), 1);
      BCL_ExampleCheck( pop.GetFitnessAverage(), 1);
      BCL_ExampleCheck( pop.GetFitnessStdDev(), 0);

      BCL_ExampleCheck( pop.GetMembers().GetSize(), 2);
      BCL_ExampleCheck( pop.GetMembersIterator().NotAtEnd(), true);

      // No uniqueness measure
      BCL_ExampleCheck( pop.HasMember( test), false);

      // With uniqueness measure
      pop.SetUniquenessMeasure( SizeTPopulationUnique());
      BCL_ExampleCheck( pop.GetUniquenessImplementation().IsDefined(), true);
      BCL_ExampleCheck( pop.HasMember( test), true);
      pop.GetUniquenessMeasure(); // call this to make sure it works
      
      // Stochastic member selection
      BCL_ExampleCheck( pop.GetRandomMember().IsDefined(), true);
      
      math::DiscreteSetSelector selector;
      BCL_ExampleCheck( pop.GetRandomMemberFromDistribution( selector).IsDefined(), true);

      // Pruning
      pop.Prune();
      BCL_ExampleCheck( pop.GetCurrentSize(), 1);

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

   }; // class ExampleOptiEvolutionPopulation

   const ExampleClass::EnumType ExampleOptiEvolutionPopulation::s_Instance
   (
     GetExamples().AddEnum( ExampleOptiEvolutionPopulation())
   );

} // namespace bcl
