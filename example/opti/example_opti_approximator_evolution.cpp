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
#include "opti/bcl_opti_approximator_evolution.h"

// includes from bcl - sorted alphabetically
#include "opti/bcl_opti_criterion_number_iterations.h"
#include "opti/bcl_opti_evolution_operation_select.h"
#include "opti/bcl_opti_evolution_population_builder.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_opti_approximator_evolution.cpp
  //! @brief this example tests the implementation of ApproximatorEvolution
  //!
  //! @author geanesar
  //! @date Nov 1, 2014
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleOptiApproximatorEvolution :
    public ExampleInterface
  {

  private:

  ////////////////////
  // Helper classes //
  ////////////////////

    //! @brief a class that builds populations by copying members from previous pops
    class PopBuilder : public opti::EvolutionPopulationBuilder< int, int>
    {
    private:

      //! single instance of this class
      static const util::SiPtr< util::ObjectInterface> s_Instance;

    public:

      //! @brief default constructor
      PopBuilder()
      {
      }

      //! @brief get the class name
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief get the class name
      const std::string &GetAlias() const
      {
        const static std::string name( "PopBuilder");
        return name;
      }

      //! @brief clone/copy function
      PopBuilder *Clone() const
      {
        return new PopBuilder( *this);
      }

      //! @brief the building method; does all the work
      void Build
      (
        opti::EvolutionPopulation< int, int> &POPULATION,
        opti::ApproximatorEvolution< int, int> &APPROXIMATOR
      )
      {
        POPULATION = opti::EvolutionPopulation< int, int>( APPROXIMATOR.GetTracker().GetCurrent()->First());
      }

      io::Serializer GetSerializer() const
      {
        return io::Serializer();
      }

    };

    //! @brief simple population scoring function
    class PopScorer :
      public opti::ApproximatorEvolution< int, int>::t_PopulationFitnessFunction
    {
    private:

      //! single instance of this class
      static const util::SiPtr< util::ObjectInterface> s_Instance;

    public:

      //! @brief constructor
      PopScorer()
      {
      }

      //! @brief copy function
      PopScorer *Clone() const
      {
        return new PopScorer( *this);
      }

      //! @brief gets the class name
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief gets the class name
      const std::string &GetAlias() const
      {
        static const std::string name( "PopScoreTest");
        return name;
      }

      //! @brief scoring function returns a constant
      int operator ()( const opti::EvolutionPopulation< int, int> &POP) const
      {
        return 1;
      }

      //! @brief serializer interface
      io::Serializer GetSerializer() const
      {
        return io::Serializer();
      }

    };

    //! @brief simple population scoring function
    class MemberScorer :
      public opti::ApproximatorEvolution< int, int>::t_MemberFitnessFunction
    {
    private:

      //! single instance of this class
      static const util::SiPtr< util::ObjectInterface> s_Instance;

    public:

      //! @brief constructor
      MemberScorer()
      {
      }

      //! @brief copy function
      MemberScorer *Clone() const
      {
        return new MemberScorer( *this);
      }

      //! @brief gets the class name
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief gets the class name
      const std::string &GetAlias() const
      {
        static const std::string name( "MemberScoreTest");
        return name;
      }

      //! @brief scoring function returns a constant
      int operator ()( const int &MEM) const
      {
        return 1;
      }

      //! @brief serializer interface
      io::Serializer GetSerializer() const
      {
        return io::Serializer();
      }

    };

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

  public:

    //! @brief Clone function
    //! @return pointer to a new ExampleOptiApproximatorEvolution
    ExampleOptiApproximatorEvolution *Clone() const
    {
      return new ExampleOptiApproximatorEvolution( *this);
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

  //////////
  // data //
  //////////

  ////////////////
  // operations //
  ////////////////

    //! @brief run routine
    //! this is performing the execution of the example
    int Run() const
    {

      // Construction and destruction
      opti::ApproximatorEvolution< int, int> approximator;

      opti::ApproximatorEvolution< int, int> approx_copy( approximator);

      // Input/output
      BCL_ExampleCheck( approximator.GetClassIdentifier(), "bcl::opti::ApproximatorEvolution<int,int>");
      BCL_ExampleCheck( approximator.GetAlias(), "Evolutionary");

      // Tracker access
      BCL_ExampleCheck( approximator.GetTrackerWithHistory().GetIteration(), 0);

      // Initial population
      opti::EvolutionPopulation< int, int> init_pop( 5);
      opti::EvolutionPopulationMember< int, int> init_mem( 1, 1);
      init_pop.AddMember( init_mem);
      init_pop.AddMember( init_mem);
      init_pop.AddMember( init_mem);
      init_pop.AddMember( init_mem);
      init_pop.AddMember( init_mem);

      // This object is declared in example_opti_evolution_operation_select.cpp
      opti::EvolutionOperationSelect< int, int> operations;
      operations.AssertRead( util::ObjectDataLabel( "(options=(Test),probabilities=(0.5))"));

      // Set/get evolution operations
      approximator.SetEvolutionOperations( operations);
      BCL_ExampleCheck( approximator.GetEvolutionOperations().GetAlias(), "OperationSelect");

      // Set/get fitness functions
      approximator.SetFitnessFunction( MemberScorer());
      BCL_ExampleCheck( approximator.GetFitnessFunction().GetAlias(), "MemberScoreTest");

      approximator.SetPopulationFitnessFunction( PopScorer());
      BCL_ExampleCheck( approximator.GetPopulationFitnessFunction().GetAlias(), "PopScoreTest");

      // Set/get inital population
      approximator.SetInitialPopulation( init_pop);
      BCL_ExampleCheck( approximator.GetTracker().GetCurrent()->First().GetCurrentSize(), 5);

      // Set the stop criterion
      opti::CriterionNumberIterations< opti::ApproximatorEvolution< int, int>::t_Population, int> criterion( 3);
      approximator.SetCriterion( criterion);

      // Set the population builder
      approximator.SetPopulationBuilder( PopBuilder());
      BCL_ExampleCheck( approximator.GetPopulationBuilder().IsDefined(), true);

      // Approximation
      approximator.Approximate();
      BCL_ExampleCheck( approximator.GetTrackerWithHistory().GetIteration(), 3);
      BCL_ExampleCheck( approximator.GetTracker().GetCurrent()->First().GetCurrentSize(), 5);

      // Do the same thing as above but with the serializer
      std::string constructor_string
      (
        "Evolutionary(population fitness=PopScoreTest,member fitness=MemberScoreTest,"
        "operations=(options=(Test),probabilities=(1.0)),builder=PopBuilder,max history=1)"
      );
      util::Implementation< opti::ApproximatorEvolution< int, int> > approx_impl;
      BCL_ExampleCheck( approx_impl.TryRead( constructor_string, util::GetLogger()), true);

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; // class ExampleOptiApproximatorEvolution

  const ExampleClass::EnumType ExampleOptiApproximatorEvolution::s_Instance
  (
    GetExamples().AddEnum( ExampleOptiApproximatorEvolution())
  );

  // Below are instance instantiations to make the serializer stuff work

  namespace opti
  {
    template<>
    const util::SiPtr< const util::ObjectInterface> ApproximatorEvolution< int, int>::s_Instance
    (
      util::Enumerated< ApproximatorEvolution< int, int> >::AddInstance
      (
        new ApproximatorEvolution< int, int>()
      )
    );
  } // namespace opti

  const util::SiPtr< util::ObjectInterface> ExampleOptiApproximatorEvolution::PopScorer::s_Instance
  (
    util::Enumerated< opti::ApproximatorEvolution< int, int>::t_PopulationFitnessFunction>::AddInstance
    (
      new ExampleOptiApproximatorEvolution::PopScorer()
    )
  );

  const util::SiPtr< util::ObjectInterface> ExampleOptiApproximatorEvolution::MemberScorer::s_Instance
  (
    util::Enumerated< opti::ApproximatorEvolution< int, int>::t_MemberFitnessFunction>::AddInstance
    (
      new ExampleOptiApproximatorEvolution::MemberScorer()
    )
  );

  const util::SiPtr< util::ObjectInterface> ExampleOptiApproximatorEvolution::PopBuilder::s_Instance
  (
    util::Enumerated< opti::EvolutionPopulationBuilder< int, int> >::AddInstance
    (
      new ExampleOptiApproximatorEvolution::PopBuilder()
    )
  );

} // namespace bcl
