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

#ifndef BCL_OPTI_PROBABILISTIC_SELECTION_H_
#define BCL_OPTI_PROBABILISTIC_SELECTION_H_

// include the namespace header
#include "bcl_opti.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_opti_approximator_evolution.h"
#include "bcl_opti_population.h"
#include "util/bcl_util_serializable_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opti
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ProbabilisticSelection
    //! @brief Interface for classes that determine when to terminate approximation
    //! @details ApproximatorModularBase derived classes will use ProbabilisticSelection derived classes to determine when
    //! to terminate their approximation. The criterion classes will be updated at each iteration and indicate when
    //! specific termination criteria are met
    //!
    //! @tparam t_MemberType is the type of the approximation argument
    //! @tparam t_FitnessType is the type of the approximation result
    //!
    //! @remarks example unnecessary
    //! @author geanesar
    //! @date May 14, 2016
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_MemberType, typename t_FitnessType>
    class ProbabilisticSelection :
      public PopulationReplacement< t_MemberType, t_FitnessType>
    {

      //! Static instance of this class
      static util::SiPtr< const util::ObjectInterface> s_Instance;

      //! the final size the population should be
      size_t m_FinalSize;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      ProbabilisticSelection
      ( 
        const size_t &FINAL_SIZE = 0
      ) :
        m_FinalSize( FINAL_SIZE)
      {
      }

      //! @brief copy the class
      //! @return a pointer to a new instance of this class
      ProbabilisticSelection *Clone() const
      {
        return new ProbabilisticSelection( *this);
      }

      //! @brief get the name for this class
      //! @return the class name
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief get the name for this class
      //! @return the class name
      const std::string &GetAlias() const
      {
        static const std::string s_name( "ProbabilisticReplacement");
        return s_name;
      }

      //! @brief selects elements from a population that will survive via tournament competition
      //! @param POPULATION the population to work on
      //! @note this will not duplicate members 
      void operator()( Population< t_MemberType, t_FitnessType> &POPULATION) const
      {
        BCL_MessageStd( "ProbabilisticReplacement operator()");

        util::SiPtrVector< const PopulationMember< t_MemberType, t_FitnessType> > members( POPULATION.GetMembers());
        size_t n_members( members.GetSize());

        if( n_members == m_FinalSize)
        {
          BCL_MessageStd( "Nothing to be done for replacement");
          return; // must use every member, so just return, don't do replacement
        }

        // fewer members than 
        BCL_Assert
        ( 
          m_FinalSize,
          "Final size for population replacement was not specified or specified as 0, cannot continue"
        );

        BCL_Assert
        ( 
          n_members > m_FinalSize, 
          "The number of members in the provided population (" + util::Format()( n_members) 
          + ") is smaller than the number required (" + util::Format()( m_FinalSize) + ")"
        );

        util::ShPtrVector< PopulationMember< t_MemberType, t_FitnessType> > kept_members;
        kept_members.AllocateMemory( m_FinalSize);

        while( !members.IsEmpty() && kept_members.GetSize() < m_FinalSize)
        {
          size_t cur_n_members( members.GetSize());
          linal::Vector< t_FitnessType> scores( cur_n_members);
          for( size_t i( 0), end_i( cur_n_members); i < end_i; ++i)
          {
            scores( i) = std::max< float>( 0.0, members( i)->GetFitness());
          }

          t_FitnessType sum( scores.Sum());
          if( sum <= 0.0)
          {
            BCL_MessageStd( "All of the scores are awful, picking a random member...");
            size_t random_member( random::GetGlobalRandom().Random< size_t>( 0, cur_n_members - 1));
            kept_members.PushBack( util::CloneToShPtr( *members( random_member)));
            members.RemoveElements( random_member, 1);
            continue;
          }

          scores /= sum;

          float test( random::GetGlobalRandom().Random< float>( 0.0, 1.0));
          size_t i( 0);
          float v( scores( i));
          for( ; v < test && i < scores.GetSize(); v += scores( i), ++i);
          if( i >= scores.GetSize())
          {
            i = scores.GetSize() - 1;
          }
          kept_members.PushBack( util::CloneToShPtr( *members( i)));
          members.RemoveElements( i, 1);
        }
        POPULATION.Reset();
        POPULATION.AddMembers( kept_members);
        BCL_MessageStd( "ProbabilisticSelection: After trimming, population now has " + util::Format()( POPULATION.GetSize()) + " members");
      }

      //! @brief hook that derived classes can override, called after TryRead successfully reads a serializer containing only initializer info (no data variables)
      //! @param SERIALIZER the serializer object with initialization information
      //! @param ERR_STREAM stream to write out errors to
      //! @return true, unless there were new errors
      bool ReadInitializerSuccessHook( const util::ObjectDataLabel &SERIALIZER, std::ostream &ERR_STREAM)
      {
        return m_FinalSize > 0;
      }

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const
      {
        io::Serializer member_data;
        member_data.SetClassDescription
        ( 
          "Determines which members of a population will be kept by running a tournament-style competition"
          " repeatedly until a particular final size is reached (i.e. a subset is sampled, and the fittest member is kept)"
        );
        member_data.AddOptionalInitializer
        (
          "final size",
          "the final size the population should be.  if greater than a population's size, or set equal to 0, nothing is done",
          io::Serialization::GetAgent( &m_FinalSize)
        );
        return member_data;
      }

    private:

    ////////////////
    // operations //
    ////////////////

    }; // template class ProbabilisticSelection< t_MemberType, t_FitnessType>
  } // namespace opti
} // namespace bcl

#endif // BCL_OPTI_PROBABILISTIC_SELECTION_H_
