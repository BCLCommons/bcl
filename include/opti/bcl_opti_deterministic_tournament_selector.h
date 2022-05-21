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

#ifndef BCL_OPTI_DETERMINISTIC_TOURNAMENT_SELECTOR_H_
#define BCL_OPTI_DETERMINISTIC_TOURNAMENT_SELECTOR_H_

// include the namespace header
#include "bcl_opti.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_opti_population.h"
#include "bcl_opti_population_member_selector.h"
#include "util/bcl_util_serializable_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opti
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class DeterministicTournamentSelector
    //! @brief Interface for classes that determine when to terminate approximation
    //! @details ApproximatorModularBase derived classes will use DeterministicTournamentSelector derived classes to determine when
    //! to terminate their approximation. The criterion classes will be updated at each iteration and indicate when
    //! specific termination criteria are met
    //!
    //! @tparam t_MemberType is the type of the approximation argument
    //! @tparam t_FitnessType is the type of the approximation result
    //!
    //! @remarks example unnecessary
    //! @author geanesar
    //! @date Sep 16, 2015
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_MemberType, typename t_FitnessType>
    class DeterministicTournamentSelector :
      public PopulationMemberSelector< t_MemberType, t_FitnessType> 
    {

      //! Static instance of this class
      static util::SiPtr< const util::ObjectInterface> s_Instance;

      //! the number of members that should be picked for each tournament
      size_t m_SampleSize;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief constructor
      DeterministicTournamentSelector( const size_t &SAMPLE = 1) :
        m_SampleSize( SAMPLE)
      {
      }

      //! @brief copy the class
      //! @return a pointer to a new instance of this class
      DeterministicTournamentSelector *Clone() const
      {
        return new DeterministicTournamentSelector( *this);
      }

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name if this function is overwritten
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief get the name of the object when used in a dynamic context
      //! @return the name of the object when used in a dynamic context
      const std::string &GetAlias() const
      {
        static const std::string s_name( "DeterministicTournament");
        return s_name;
      }

    ///////////////
    // operators //
    ///////////////

      //! @brief competes a subset of the population against each other and takes the fittest one
      //! @note if the sample size is 1, this randomly selects a member from the population
      //! @return a SiPtr to the chosen member 
      util::SiPtr< const PopulationMember< t_MemberType, t_FitnessType> > operator()
      (
        const Population< t_MemberType, t_FitnessType> &POPULATION
      ) const
      {

        util::SiPtr< const PopulationMember< t_MemberType, t_FitnessType> > winner;
        util::SiPtrVector< const PopulationMember< t_MemberType, t_FitnessType> > members( POPULATION.GetMembers());
        if( members.IsEmpty())
        {
          BCL_MessageStd( "Population was empty; cannot select any members via tournament");
          return winner;
        }

        // Since the tournament will ultimately take the highest fitness
        size_t n_needed( m_SampleSize);
        size_t n_available( members.GetSize());
        for( size_t i( 0), end_i( members.GetSize()); i < end_i && n_needed < n_available; --n_available, ++i)
        {

          float prob( float( n_needed) / float( n_available));
          float test( random::GetGlobalRandom().Random< float>( 0.0, 1.0));
          if( test < prob)
          {
            --n_needed;
            if( !winner.IsDefined() || winner->GetFitness() < members( i)->GetFitness())
            {
              winner = members( i);
            }
          }
        }
        return winner;
      }

      //! @brief hook that derived classes can override, called after TryRead successfully reads a serializer containing only initializer info (no data variables)
      //! @param SERIALIZER the serializer object with initialization information
      //! @param ERR_STREAM stream to write out errors to
      //! @return true, unless there were new errors
      bool ReadInitializerSuccessHook( const util::ObjectDataLabel &SERIALIZER, std::ostream &ERR_STREAM)
      {
        if( !m_SampleSize)
        {
          BCL_MessageStd( "Tournament cannot be used with a zero tournament size");
        }
        return m_SampleSize > 0;
      }

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const
      {
        io::Serializer member_data;
        member_data.SetClassDescription
        ( 
          "Selects population members in a deterministic tournament, i.e. "
          "selects a subset of the population and selects the fittest member"
        );
        member_data.AddInitializer
        (
          "tournament size",
          "the number of members that should be compared during each tournament (1 is equivalent to random selection)",
          io::Serialization::GetAgent( &m_SampleSize),
          "1"
        );
        return member_data;
      }

    }; // template class DeterministicTournamentSelector< t_MemberType, t_FitnessType>
  } // namespace opti
} // namespace bcl

#endif // BCL_OPTI_DETERMINISTIC_TOURNAMENT_SELECTOR_H_
