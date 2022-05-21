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
#include "opti/bcl_opti_evolution_operation_select.h"

// includes from bcl - sorted alphabetically
#include "opti/bcl_opti_evolution_population.hpp"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_opti_evolution_operation_select.cpp
  //! @brief this example tests the implementation of EvolutionOperationSelect
  //!
  //! @author geanesar
  //! @date Nov 1, 2014
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleOptiEvolutionOperationSelect :
    public ExampleInterface
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

  public:

    //! @brief Clone function
    //! @return pointer to a new ExampleOptiEvolutionOperationSelect
    ExampleOptiEvolutionOperationSelect *Clone() const
    {
      return new ExampleOptiEvolutionOperationSelect( *this);
    }

  //////////
  // data //
  //////////

    //! Typedefs
    typedef opti::ApproximatorEvolution< int, int>::EvolutionOperationSerializable EvolutionOperationSerializable;

    //! @brief a simple evolution operation
    class EvolutionOperation :
      public opti::ApproximatorEvolution< int, int>::EvolutionOperationSerializable
    {
    private:

      //! static instance of this class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    public:

      //! @brief default constructor
      EvolutionOperation()
      {
      }

      //! @brief Clone function
      //! @return a pointer to a copy of this class
      EvolutionOperation *Clone() const
      {
        return new EvolutionOperation( *this);
      }

      //! @brief the name of this class
      //! @return the name of this class
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( this);
      }

      //! @brief the name of this class
      //! @return the name of this class
      const std::string &GetAlias() const
      {
        const static std::string name( "Test");
        return name;
      }

      //! @brief operator to mutate an integer by adding 1 to it
      //! @param FEED member feed to use to get the next integer
      //! @return a vector with a single integer in it
      util::ShPtrVector< opti::EvolutionPopulationMember< int, int> > operator()
      (
        const util::SiPtrVector< const opti::EvolutionPopulationMember< int, int> > &MEMBERS
      ) const
      {
        util::ShPtrVector< opti::EvolutionPopulationMember< int, int> > return_vector;
        for( size_t i( 0); i < MEMBERS.GetSize(); ++i)
        {
          return_vector.PushBack
            (
             util::ShPtr< opti::EvolutionPopulationMember< int, int> >
             (
              new opti::EvolutionPopulationMember< int, int>( MEMBERS( i)->GetMember() + 1, MEMBERS( i)->GetMember() + 1)
             )
            );
        }
        return return_vector;
      }

      //! @brief gets the input size
      //! @return input size
      size_t GetInputSize() const
      {
        return 0;
      }

      //! @brief gets information about member data for serialization
      //! @return a serializer which contains info about member data
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

      // Constructor
      opti::EvolutionOperationSelect< int, int> selector;

      // Cloning
      util::ShPtr< opti::EvolutionOperationSelect< int, int> > sp_selector( selector.Clone());
      BCL_ExampleCheck( sp_selector.IsDefined(), true);

      // Serialized construction
      std::string constructor( "(options=(Test),probabilities=(0.5))");
      BCL_ExampleCheck( sp_selector->TryRead( constructor, util::GetLogger()), true);

      // Destruction
      sp_selector.Reset();

      // Add an operation
      selector.AddImplementation( EvolutionOperation(), 0.5);
      BCL_ExampleIndirectCheck( selector.GetSize(), 1, "addition of an implementation");

      //! Get the decision node
      math::ObjectStochasticSelector< EvolutionOperationSerializable> &decision_node( selector.GetDecisionNode());

      BCL_ExampleIndirectCheck( decision_node.GetSize(), 1, "fetching the decision node works");

      const EvolutionOperationSerializable &chosen( selector.Select());
      BCL_ExampleIndirectCheck( chosen.GetAlias(), "Test", "selector chose an option");

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;
    
  }; // class ExampleOptiEvolutionOperationSelect

  const ExampleClass::EnumType ExampleOptiEvolutionOperationSelect::s_Instance
  (
    GetExamples().AddEnum( ExampleOptiEvolutionOperationSelect())
  );

  const util::SiPtr< const util::ObjectInterface> ExampleOptiEvolutionOperationSelect::EvolutionOperation::s_Instance
  (
    util::Enumerated< opti::ApproximatorEvolution< int, int>::EvolutionOperationSerializable>::AddInstance
    (
      new ExampleOptiEvolutionOperationSelect::EvolutionOperation()
    )
  );

} // namespace bcl
