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

#ifndef BCL_OPTI_APPROXIMATOR_EVOLUTION_H_
#define BCL_OPTI_APPROXIMATOR_EVOLUTION_H_

// include the namespace header
#include "bcl_opti.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_opti_approximator_modular_base.h"
#include "bcl_opti_evolution_operation_select.h"
#include "bcl_opti_evolution_population.h"
#include "bcl_opti_evolution_population_builder.h"
#include "bcl_opti_evolution_population_member.h"
#include "bcl_opti_tracker_with_history.h"
#include "math/bcl_math_function_interface_serializable.h"
#include "util/bcl_util_implementation.h"
#include "util/bcl_util_message.h"
#include "util/bcl_util_serializable_interface.h"
#include "util/bcl_util_sh_ptr_vector.h"
#include "util/bcl_util_undefined.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opti
  {

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ApproximatorEvolution
    //! @brief Implementation of an approximator for an evolutionary algorithm
    //!
    //! @tparam t_MemberType is the class that will be modified/evolved
    //! @tparam t_FitnessType is the type of fitness the class will use
    //!
    //! @see http://en.wikipedia.org/wiki/Evolutionary_algorithm
    //!
    //! @see @link example_opti_approximator_evolution.cpp @endlink
    //! @author geanesar
    //! @date Oct 29, 2014
    //!
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_MemberType, typename t_FitnessType>
    class ApproximatorEvolution :
      public ApproximatorModularBase< EvolutionPopulation< t_MemberType, t_FitnessType>, t_FitnessType>,
      public util::SerializableInterface
    {
    
    public:

      //! A single instance of this class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////
    // typedefs //
    //////////////

      //! Typedef for populations
      typedef EvolutionPopulation< t_MemberType, t_FitnessType> t_Population;

      //! Population members
      typedef EvolutionPopulationMember< t_MemberType, t_FitnessType> t_PopulationMember;

      //! Population scoring function
      typedef math::FunctionInterfaceSerializable< t_Population, t_FitnessType> t_PopulationFitnessFunction;

      //! Member scoring function
      typedef math::FunctionInterfaceSerializable< t_MemberType, t_FitnessType> t_MemberFitnessFunction;

    //////////
    // data //
    //////////

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //!
      //! @class EvolutionOperationSerializable
      //! @brief Simple aliasing class for cleaning up command line help (effective typedef)
      //! @details Classes that derive from this one should output a set of EvolutionPopulationMembers that have their
      //! @details histories set so that the parents can be identified later on
      //!
      //! @remark example unnecessary
      //! @author geanesar
      //! @date Oct 29, 2014
      //!
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      class EvolutionOperationSerializable :
        public math::FunctionInterfaceSerializable
        <
          util::SiPtrVector< const t_PopulationMember>,
          util::ShPtrVector< t_PopulationMember>
        > 
      {
      public:

        //! @brief clone function for copying the class
        //! @return a pointer to a copy of this class
        virtual EvolutionOperationSerializable *Clone() const = 0;

        //! @brief gets the number of members needed as input
        //! @return the number needed
        virtual size_t GetInputSize() const = 0;
      };
        
    private:

      //! The tracker; stored to reduce casting overhead
      util::SiPtr< TrackerWithHistory< t_Population, t_FitnessType> > m_TrackerWithHistory;

      //! Fitness function which scores an entire population
      util::Implementation< math::FunctionInterfaceSerializable< t_Population, t_FitnessType> > m_PopulationFitnessFunction;

      //! Function to score member fitness
      util::Implementation< math::FunctionInterfaceSerializable< t_MemberType, t_FitnessType> > m_MemberFitnessFunction;

      //! Available evolutionary operations
      EvolutionOperationSelect< t_MemberType, t_FitnessType> m_EvolutionOperations;

      //! The population builder
      util::Implementation< EvolutionPopulationBuilder< t_MemberType, t_FitnessType> > m_PopulationBuilder;

      //! histories list
      storage::List< util::ObjectDataLabel> m_HistoriesList;

      //! The number of populations to keep track of
      size_t m_NumberPopsToTrack;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      ApproximatorEvolution() :
        ApproximatorModularBase< t_Population, t_FitnessType>
        (
          TrackerWithHistory< t_Population, t_FitnessType>( e_LargerIsBetter)
        ),
        m_TrackerWithHistory( &( this->GetTracker()))
      {
      }

      //! @brief copy constructor, necessary to copy tracker pointers and implementations
      ApproximatorEvolution( const ApproximatorEvolution< t_MemberType, t_FitnessType> &OTHER) :
        ApproximatorModularBase< t_Population, t_FitnessType>( OTHER),
        m_TrackerWithHistory( &( this->GetTracker())),
        m_PopulationFitnessFunction( OTHER.m_PopulationFitnessFunction),
        m_MemberFitnessFunction( OTHER.m_MemberFitnessFunction),
        m_EvolutionOperations( OTHER.m_EvolutionOperations),
        m_NumberPopsToTrack( OTHER.m_NumberPopsToTrack)
      {
        this->GetTrackerWithHistory().SetMaxHistorySize( m_NumberPopsToTrack);
      }

      //! @brief Clone function
      //! @return pointer to a new ApproximatorEvolution< t_MemberType>
      ApproximatorEvolution *Clone() const
      {
        return new ApproximatorEvolution( *this);
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

      //! @brief returns the class name
      //! @return the class name when used in a dynamic context
      const std::string &GetAlias() const
      {
        static std::string s_alias( "Evolutionary");
        return s_alias;
      }

      //! @brief get a tracker with history from this class
      //! @return the tracker with history (non-const reference)
      TrackerWithHistory< t_Population, t_FitnessType> &GetTrackerWithHistory()
      {
        return *m_TrackerWithHistory;
      }

      //! @brief get a tracker with history from this class
      //! @return the tracker with history (const reference)
      const TrackerWithHistory< t_Population, t_FitnessType> &GetTrackerWithHistory() const
      {
        return *m_TrackerWithHistory;
      }

      //! @brief sets the initial population and tracks it
      //! @param POPULATION the population to use
      //! @param RESCORE whether to rescore the EvolutionPopulationMemberFilterInterfacemembers of POPULATION
      void SetInitialPopulation
      (
        const t_Population &POPULATION,
        const bool &RESCORE = true
      )
      {
        // Reset the tracker
        this->m_HistoriesList.Reset();
        this->GetTrackerWithHistory().Reset();
        this->GetTrackerWithHistory().SetPhase( e_Start);
        
        t_FitnessType initial_fitness( m_PopulationFitnessFunction->operator()( POPULATION));

        //! Copy the initial population
        util::ShPtr< storage::Pair< t_Population, t_FitnessType> > sp_initial_population
        (
          new storage::Pair< t_Population, t_FitnessType>
          (
            t_Population( POPULATION),
            initial_fitness
          )
        );

        this->Track( sp_initial_population);
      }

      //! @brief return the histories list
      //! @return a list containing datalabels that represent histories
      const storage::List< util::ObjectDataLabel> &GetHistoryList() const
      {
        return m_HistoriesList;
      }

      //! @brief get the fitness function used to score members in the algorithm
      //! @return the member fitness function
      const t_MemberFitnessFunction &GetFitnessFunction() const
      {
        return *m_MemberFitnessFunction;
      }

      //! @brief set the member fitness function
      //! @param FUNCTION the fitness function to use
      void SetFitnessFunction
      (
        const t_MemberFitnessFunction &FUNCTION
      )
      {
        m_MemberFitnessFunction = FUNCTION.GetLabel();
      }

      //! @brief get the fitness function used to score populations
      //! @return the population fitness function
      const t_PopulationFitnessFunction &GetPopulationFitnessFunction() const
      {
        return *m_PopulationFitnessFunction;
      }

      //! @brief set the population fitness function
      //! @param FUNCTION the function to use
      void SetPopulationFitnessFunction
      (
        const math::FunctionInterfaceSerializable< t_Population, t_FitnessType> &FUNCTION
      )
      {
        m_PopulationFitnessFunction = FUNCTION.GetLabel();
      }

      //! @brief set the population builder
      //! @param BUILDER the builder to use
      void SetPopulationBuilder( const EvolutionPopulationBuilder< t_MemberType, t_FitnessType> &BUILDER)
      {
        m_PopulationBuilder = util::Implementation< EvolutionPopulationBuilder< t_MemberType, t_FitnessType> >( BUILDER);
      }

      //! @brief get the population builder
      //! @return a reference to the population builder implementation
      const util::Implementation< EvolutionPopulationBuilder< t_MemberType, t_FitnessType> >
        &GetPopulationBuilder() const
      {
        return m_PopulationBuilder;
      }

      //! @brief get the evolution operations that the algorithm will use
      //! @return the evolution operations
      const EvolutionOperationSelect< t_MemberType, t_FitnessType> &GetEvolutionOperations() const
      {
        return m_EvolutionOperations;
      }

      //! @brief set the evolution operations to use
      //! @param OPERATIONS the operations to use
      void SetEvolutionOperations
      (
        const EvolutionOperationSelect< t_MemberType, t_FitnessType> &OPERATIONS
      )
      {
        m_EvolutionOperations = OPERATIONS;
      }

    private:

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief sets the name on the ObjectDataLabel used for the member's history
      //! @param MEMBER the member to modify
      //! @param ITERATION the generation number
      //! @param MEMBER_NUMBER the number of the member in the population (arbitrary number, not its rank)
      void AddIdentifierAndFitnessToHistory
      (
        t_PopulationMember &MEMBER,
        const size_t &ITERATION,
        const size_t &MEMBER_NUMBER
      ) const
      {

        // the ID/name of this member, unique amongst all members made from this approximator
        std::string identifier( "P" + util::Format()( ITERATION) + "M" + util::Format()( MEMBER_NUMBER));

        // If the member does not have a "fitness" field in its history already, go ahead and set one
        storage::Vector< util::ObjectDataLabel> new_history( MEMBER.GetHistory().GetArguments());
        if( MEMBER.GetHistory().FindName( "fitness") == MEMBER.GetHistory().End())
        {
          new_history.PushBack( util::ObjectDataLabel( "fitness", util::Format()( MEMBER.GetFitness())));
        }

        // Strip the original identifier, replace it with the correct one of the format
        // PXMY=(fitness=<fitness>,<remaining stuff>)
        MEMBER.SetHistory( util::ObjectDataLabel( identifier, "", new_history));
      }

    protected:

      //! @brief conducts the next approximation step and stores the approximation
      void Next()
      {

        // Make a new pair that we can track later on, so that we don't have to copy a population later
        util::ShPtr< storage::Pair< t_Population, t_FitnessType> > sp_next_pop
        (
          new storage::Pair< t_Population, t_FitnessType>
          (
            t_Population( this->GetTracker().GetCurrent()->First().GetNormalSize()),
            util::GetUndefined< t_FitnessType>()
          )
        );

        // reference for convenience
        t_Population &next_population
        (
          sp_next_pop->First()
        );

        // Set uniqueness measure
        const util::Implementation< EvolutionMemberUniqueInterface< t_MemberType, t_FitnessType> > &uniqueness
        (
          this->GetTracker().GetCurrent()->First().GetUniquenessImplementation()
        );

        if( uniqueness.IsDefined())
        {
          next_population.SetUniquenessMeasure( *uniqueness);
        }

        // Build the population using specified method 
        BCL_Assert
        (
          m_PopulationBuilder.IsDefined(),
          "Population builder is not defined"
        );
        m_PopulationBuilder->Build( next_population, *this);

        // Score the population
        if( m_PopulationFitnessFunction.IsDefined())
        {
          sp_next_pop->Second() = ( *m_PopulationFitnessFunction)( next_population);
        }

        // Set the member histories appropriately
        iterate::Generic< t_PopulationMember> itr_member
        (
          next_population.GetMembersIterator()
        );

        for( size_t mem_no( 0); itr_member.NotAtEnd(); ++itr_member, ++mem_no)
        {
          AddIdentifierAndFitnessToHistory( *itr_member, this->GetTracker().GetIteration(), mem_no);
          m_HistoriesList.PushBack
          (
            itr_member->GetHistory()
          );
        }

        this->Track
        (
          sp_next_pop
        );
        BCL_MessageStd( "Finished population " + util::Format()( this->GetTracker().GetIteration()));
      }

      //! @brief evaluates whether the approximation can continue
      //! @return true, if the approximation can continue - otherwise false
      bool CanContinue() const
      {
        return true;
      }

    //////////////////////
    // input and output //
    //////////////////////

    public:

      //! @brief called after TryRead successfully reads a serializer containing only initializer info (no data variables)
      //! @param SERIALIZER the serializer object with initialization information
      //! @param ERR_STREAM stream to write out errors to
      //! @return true, unless there were new errors
      bool ReadInitializerSuccessHook( const util::ObjectDataLabel &SERIALIZER, std::ostream &ERR_STREAM)
      {
        this->GetTrackerWithHistory().SetMaxHistorySize( m_NumberPopsToTrack);
        return m_PopulationFitnessFunction.IsDefined()
          && m_MemberFitnessFunction.IsDefined();
      }

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const
      {
        io::Serializer member_data;

        member_data.AddInitializer
        (
          "population fitness",
          "fitness function for scoring members",
          io::Serialization::GetAgent( &m_PopulationFitnessFunction),
          ""
        );

        member_data.AddInitializer
        (
          "member fitness",
          "fitness function for scoring members",
          io::Serialization::GetAgent( &m_MemberFitnessFunction),
          ""
        );

        member_data.AddInitializer
        (
          "operations",
          "operations that will be used to mutate members",
          io::Serialization::GetAgent( &m_EvolutionOperations),
          ""
        );

        member_data.AddInitializer
        (
          "builder",
          "how populations will be built",
          io::Serialization::GetAgent( &m_PopulationBuilder),
          "Simple"
        );

        member_data.AddInitializer
        (
          "max history",
          "the maximum number of populations to keep track of.  this should be specified if the population builder requires "
          "more than one population.  set this to 0 to track every population (this may consume lots of memory)",
          io::Serialization::GetAgent( &m_NumberPopsToTrack),
          "1"
        );

        member_data.SetClassDescription( "An implementation of an evolutionary algorithm for optimization problems");

        return member_data;
      }
        
      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM)
      {
        // end
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        // end
        return OSTREAM;
      }

    private:

      //! @brief an operation nobody should ever use
      ApproximatorEvolution &operator =( const ApproximatorEvolution &);

    }; // template class ApproximatorEvolution< t_MemberType, t_FitnessType>

  } // namespace opti
} // namespace bcl

#endif // BCL_OPTI_APPROXIMATOR_EVOLUTION_H_
