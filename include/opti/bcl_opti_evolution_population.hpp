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

#ifndef BCL_OPTI_EVOLUTION_POPULATION_HPP_
#define BCL_OPTI_EVOLUTION_POPULATION_HPP_

// include the namespace header
#include "bcl_opti.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_opti_evolution_population.h"
#include "bcl_opti_evolution_population_member.h"
#include "linal/bcl_linal_vector.h"
#include "math/bcl_math_discrete_set_selector.h"
#include "math/bcl_math_running_average_sd.h"
#include "util/bcl_util_implementation.h"
#include "util/bcl_util_ptr_interface.h"
#include "util/bcl_util_sh_ptr_vector.h"
#include "util/bcl_util_si_ptr_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opti
  {

    //! @brief sorts members by fitness score
    template< typename t_MemberType, typename t_FitnessType>
    void EvolutionPopulation< t_MemberType, t_FitnessType>::SortMembers()
    {
      m_FitnessVectorStale = true;
      m_Members.Sort
      (
        EvolutionMemberCompareLess< t_MemberType, t_FitnessType>()
      );
    }

    //! @brief Recalculates the fitness statistics (average, std dev)
    template< typename t_MemberType, typename t_FitnessType>
    void EvolutionPopulation< t_MemberType, t_FitnessType>::RecalculateFitnessStatistics()
    {
      m_FitnessAvg.Reset();
      for( size_t pos( 0); pos < m_Members.GetSize(); ++pos)
      {
        m_FitnessAvg += m_Members( pos)->GetFitness();
      }
    }

    //! @brief regenerates m_FitnessVector from current members
    template< typename t_MemberType, typename t_FitnessType>
    void EvolutionPopulation< t_MemberType, t_FitnessType>::RegenerateFitnessVector() const
    {
      linal::Vector< t_FitnessType> new_fitness_vector( m_Members.GetSize());

      for( size_t pos( 0); pos < m_Members.GetSize(); ++pos)
      {
        new_fitness_vector( pos) = m_Members( pos)->GetFitness();
      }

      m_FitnessVector = new_fitness_vector;
      m_FitnessVectorStale = false;
    }

    //! @brief Constructor with unspecified normal size
    template< typename t_MemberType, typename t_FitnessType>
    EvolutionPopulation< t_MemberType, t_FitnessType>::EvolutionPopulation() :
      m_NormalPopulationSize( 0),
      m_Members(),
      m_MemberUnique()
    {
    }

    //! @brief constructor with population size
    template< typename t_MemberType, typename t_FitnessType>
    EvolutionPopulation< t_MemberType, t_FitnessType>::EvolutionPopulation( const size_t &POPULATION_SIZE) :
      m_NormalPopulationSize( POPULATION_SIZE),
      m_Members(),
      m_MemberUnique()
    {
      m_Members.AllocateMemory( POPULATION_SIZE);
    }

    //! @brief Constructor with population size and initial memberes
    //! @param POPULATION_SIZE the desired size of the population
    //! @param INITIAL_MEMBERS the members to add immediately
    template< typename t_MemberType, typename t_FitnessType>
    EvolutionPopulation< t_MemberType, t_FitnessType>::EvolutionPopulation
    (
      const size_t &POPULATION_SIZE,
      util::ShPtrVector< EvolutionPopulationMember< t_MemberType, t_FitnessType> > &INITIAL_MEMBERS
    ) :
      m_NormalPopulationSize( POPULATION_SIZE),
      m_Members( INITIAL_MEMBERS),
      m_MemberUnique()
    {
      RecalculateFitnessStatistics();
    }

    //! @brief copy constructor, does a deep copy of the members
    template< typename t_MemberType, typename t_FitnessType>
    EvolutionPopulation< t_MemberType, t_FitnessType>::EvolutionPopulation
    (
      const EvolutionPopulation< t_MemberType, t_FitnessType> &OTHER
    ) :
      m_NormalPopulationSize( OTHER.m_NormalPopulationSize),
      m_MemberUnique( OTHER.m_MemberUnique)
    {
      size_t final_size( OTHER.m_Members.GetSize());
      m_Members.AllocateMemory( final_size);
      for( size_t i( 0); i < final_size; ++i)
      {
        m_Members.PushBack( util::CloneToShPtr( *OTHER.m_Members( i)));
      }
      RegenerateFitnessVector();
      m_FitnessVectorStale = false;
      RecalculateFitnessStatistics();
    }

    //! @brief Clone function
    //! @return A copy of this class
    template< typename t_MemberType, typename t_FitnessType> EvolutionPopulation< t_MemberType, t_FitnessType> *
      EvolutionPopulation< t_MemberType, t_FitnessType>::Clone() const
    {
      return new EvolutionPopulation( *this);
    }

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    template< typename t_MemberType, typename t_FitnessType>
    const std::string &EvolutionPopulation< t_MemberType, t_FitnessType>::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief Get the members from the population
    //! @return Returns a const reference to all members
    template< typename t_MemberType, typename t_FitnessType>
    util::SiPtrVector< const EvolutionPopulationMember< t_MemberType, t_FitnessType> >
    EvolutionPopulation< t_MemberType, t_FitnessType>::GetMembers() const
    {
      return util::SiPtrVector< const EvolutionPopulationMember< t_MemberType, t_FitnessType> >
      (
        m_Members
      );
    }

    //! @brief gets an iterator for members
    //! @return a generic iterator to members of the populations
    template< typename t_MemberType, typename t_FitnessType>
    iterate::Generic< const EvolutionPopulationMember< t_MemberType, t_FitnessType> >
    EvolutionPopulation< t_MemberType, t_FitnessType>::GetMembersIterator() const
    {
      return iterate::Generic< const EvolutionPopulationMember< t_MemberType, t_FitnessType> >
      (
        m_Members.Begin(), m_Members.End()
      );
    }

    //! @brief gets an iterator for members
    //! @return a generic iterator to members of the populations
    template< typename t_MemberType, typename t_FitnessType>
    iterate::Generic< EvolutionPopulationMember< t_MemberType, t_FitnessType> >
    EvolutionPopulation< t_MemberType, t_FitnessType>::GetMembersIterator()
    {
      return iterate::Generic< EvolutionPopulationMember< t_MemberType, t_FitnessType> >
      (
        m_Members.Begin(), m_Members.End()
      );
    }

    //! @brief Return the top percentage of fittest members
    //! @param PERCENTAGE the percentage of members to return, starting with the fittest
    //! @return SiPtrVector containing pairs of members and their fitness scores
    template< typename t_MemberType, typename t_FitnessType>
    util::SiPtrVector< const EvolutionPopulationMember< t_MemberType, t_FitnessType> >
    EvolutionPopulation< t_MemberType, t_FitnessType>::GetFittestPercent
    (
      const float &PERCENTAGE
    ) const
    {
      if( PERCENTAGE < 0.0 || PERCENTAGE > 1.0)
      {
        return util::SiPtrVector< const EvolutionPopulationMember< t_MemberType, t_FitnessType> >();
      }
      size_t members_to_return( m_Members.GetSize() * PERCENTAGE);
      return GetFittestMembers( members_to_return);
    }

    //! @brief Get a number of the fittest members of a population
    //! @param NUMBER the number of members to return, starting with the fittest
    //! @return SiPtrVector to fittest members
    template< typename t_MemberType, typename t_FitnessType>
    util::SiPtrVector< const EvolutionPopulationMember< t_MemberType, t_FitnessType> >
    EvolutionPopulation< t_MemberType, t_FitnessType>::GetFittestMembers
    (
      const size_t &NUMBER
    ) const
    {
      // Sort members to ease this process
      SortMembers();

      // Add fittest members to the population
      size_t num_members( NUMBER < m_Members.GetSize() ? NUMBER : m_Members.GetSize());
      util::SiPtrVector< const EvolutionPopulationMember< t_MemberType, t_FitnessType> > return_vector( num_members);
      size_t num_added( 0);
      for( size_t pos( 0); pos < num_members; ++pos)
      {
        return_vector( num_added) = m_Members( pos);
      }
      return return_vector;
    }

    //! @brief get members with a fitness over a cutoff
    //! @param FITNESS minimum fitness
    //! @return SiPtrVector
    template< typename t_MemberType, typename t_FitnessType>
    util::SiPtrVector< const EvolutionPopulationMember< t_MemberType, t_FitnessType> >
    EvolutionPopulation< t_MemberType, t_FitnessType>::GetMembersWithMinimumFitness
    (
      const t_FitnessType &FITNESS
    ) const
    {
      util::SiPtrVector< const EvolutionPopulationMember< t_MemberType, t_FitnessType> > return_vector;

      // Add members to the vector
      for( size_t pos( 0); pos < m_Members.GetSize(); ++pos)
      {
        if( m_Members( pos)->GetFitness() >= FITNESS)
        {
          return_vector.PushBack( m_Members( pos));
        }
      }
      return return_vector;
    }

    //! @brief return a member using a probability distribution
    //! @param SELECTOR the discrete set selector (i.e. probability distribution function) to use
    //! @return SiPtr to the chosen member
    template< typename t_MemberType, typename t_FitnessType>
    util::SiPtr< const EvolutionPopulationMember< t_MemberType, t_FitnessType> >
    EvolutionPopulation< t_MemberType, t_FitnessType>::GetRandomMemberFromDistribution
    (
      const math::DiscreteSetSelector &SELECTOR
    ) const
    {
      // Copy the selector (since we need to add things to it
      math::DiscreteSetSelector selector( SELECTOR);
      selector.Reset();
      if( m_FitnessVectorStale)
      {
        RegenerateFitnessVector();
      }
      selector.Prepare( linal::Vector< double>( m_FitnessVector.Begin(), m_FitnessVector.End()));
      size_t pos( selector());
      //BCL_MessageStd( "Returning member " + util::Format()( pos));
      return util::SiPtr< const EvolutionPopulationMember< t_MemberType, t_FitnessType> >( m_Members( pos));
    }

    //! @brief Get a random member of the population (uniform probability)
    //! @return a pointer to the chosen member
    template< typename t_MemberType, typename t_FitnessType>
    util::SiPtr< const EvolutionPopulationMember< t_MemberType, t_FitnessType> >
    EvolutionPopulation< t_MemberType, t_FitnessType>::GetRandomMember() const
    {
      return GetRandomMemberFromDistribution( math::DiscreteSetSelector());
    }

    //! @brief Get the current size of the population
    //! @return The number of members currently held in the population
    template< typename t_MemberType, typename t_FitnessType>
    size_t EvolutionPopulation< t_MemberType, t_FitnessType>::GetCurrentSize() const
    {
      return m_Members.GetSize();
    }

    //! @brief Get the normal size of the population, i.e. the size of a full population
    //! @return the size of a full population
    template< typename t_MemberType, typename t_FitnessType>
    size_t EvolutionPopulation< t_MemberType, t_FitnessType>::GetNormalSize() const
    {
      return m_NormalPopulationSize;
    }

    //! @brief set the normal size of the population
    //! @param SIZE the size
    template< typename t_MemberType, typename t_FitnessType>
    void EvolutionPopulation< t_MemberType, t_FitnessType>::SetNormalSize( const size_t &SIZE)
    {
      m_NormalPopulationSize = SIZE;
    }

    //! @brief Clears the population (deletes members)
    template< typename t_MemberType, typename t_FitnessType>
    void EvolutionPopulation< t_MemberType, t_FitnessType>::Reset()
    {
      m_Members.Reset();
      m_FitnessAvg.Reset();
      m_FitnessVector = linal::Vector< float>();
      m_FitnessVectorStale = true;
    }

    //! @brief Check if the population has a particular member, then add it
    //! @param MEMBER the member to add
    //! @return true if the member was added, false otherwise
    template< typename t_MemberType, typename t_FitnessType>
    bool EvolutionPopulation< t_MemberType, t_FitnessType>::HasMember
    (
      const EvolutionPopulationMember< t_MemberType, t_FitnessType> &MEMBER
    ) const
    {
      return m_MemberUnique.IsDefined() ? m_MemberUnique->Contains( MEMBER, *this) : false;
    }

    //! @brief Gets the lowest fitness present in the population
    //! @return The value of the lowest fitness in the population, or undefined if the population is empty
    template< typename t_MemberType, typename t_FitnessType>
    t_FitnessType EvolutionPopulation< t_MemberType, t_FitnessType>::GetLowestFitness() const
    {
      if( m_Members.GetSize() == 0)
      {
        return util::GetUndefined< t_FitnessType>();
      }

      t_FitnessType lowest( m_Members( 0)->GetFitness());
      for( size_t pos( 0); pos < m_Members.GetSize(); ++pos)
      {
        if( m_Members( pos)->GetFitness() < lowest)
        {
          lowest = m_Members( pos)->GetFitness();
        }
      }
      return lowest;
    }

    //! @brief Gets the highest fitness in the population
    //! @return The value of the highest fitness in the population, or undefined if the population is empty
    template< typename t_MemberType, typename t_FitnessType>
    t_FitnessType EvolutionPopulation< t_MemberType, t_FitnessType>::GetHighestFitness() const
    {
      if( m_Members.GetSize() == 0)
      {
        return util::GetUndefined< t_FitnessType>();
      }

      t_FitnessType highest( m_Members( 0)->GetFitness());
      for( size_t pos( 0); pos < m_Members.GetSize(); ++pos)
      {
        if( m_Members( pos)->GetFitness() > highest)
        {
          highest = m_Members( pos)->GetFitness();
        }
      }
      return highest;
    }

    //! @brief Gets the average fitness in the population
    //! @return The value of the average fitness of the population
    template< typename t_MemberType, typename t_FitnessType>
    t_FitnessType EvolutionPopulation< t_MemberType, t_FitnessType>::GetFitnessAverage() const
    {
      return m_FitnessAvg.GetAverage();
    }

    //! @brief Gets the standard deviation of the fitnesses
    //! @return the fitness standard deviation
    template< typename t_MemberType, typename t_FitnessType>
    t_FitnessType EvolutionPopulation< t_MemberType, t_FitnessType>::GetFitnessStdDev() const
    {
      return m_FitnessAvg.GetStandardDeviation();
    }

    //! @brief Adds a member to the population
    //! @param MEMBER the member to add
    //! @param FITNESS_SCORE the fitness of the new member
    //! @return true if the member was added successfully
    template< typename t_MemberType, typename t_FitnessType>
    bool EvolutionPopulation< t_MemberType, t_FitnessType>::AddMember
    (
      const t_MemberType &MEMBER,
      const t_FitnessType &FITNESS_SCORE
    )
    {
      return AddMember
        (
          util::ShPtr< EvolutionPopulationMember< t_MemberType, t_FitnessType> >
          (
            new EvolutionPopulationMember< t_MemberType, t_FitnessType>( MEMBER, FITNESS_SCORE)
          )
        );
    }

    //! @brief add a member to the population, then re-sort the members
    //! @param NEW_MEMBER the member to add
    template< typename t_MemberType, typename t_FitnessType>
    bool EvolutionPopulation< t_MemberType, t_FitnessType>::AddMember
    (
      const EvolutionPopulationMember< t_MemberType, t_FitnessType> &NEW_MEMBER
    )
    {
      return AddMember( util::CloneToShPtr( NEW_MEMBER));
    }

    //! @brief add a member to the population
    //! @param NEW_MEMBER the member to add
    template< typename t_MemberType, typename t_FitnessType>
    bool EvolutionPopulation< t_MemberType, t_FitnessType>::AddMember
    (
      const util::ShPtr< EvolutionPopulationMember< t_MemberType, t_FitnessType> > &NEW_MEMBER
    )
    {

      if( !this->HasMember( *NEW_MEMBER))
      {
        m_FitnessVectorStale = true;
        m_Members.PushBack( NEW_MEMBER);
        m_FitnessAvg += NEW_MEMBER->GetFitness();

        //SortMembers();
        return true;
      }
      return false;
    }

    //! @brief add members to the population
    //! @param MEMBERS the members to add
    template< typename t_MemberType, typename t_FitnessType>
    size_t EvolutionPopulation< t_MemberType, t_FitnessType>::AddMembers( const util::ShPtrVector< EvolutionPopulationMember< t_MemberType, t_FitnessType> > &MEMBERS)
    {
      size_t num_added( 0);
      m_Members.AllocateMemory( MEMBERS.GetSize() + m_Members.GetSize());
      for( size_t i( 0), end( MEMBERS.GetSize()); i < end; ++i)
      {
        if( AddMember( *MEMBERS( i)))
        {
          ++num_added;
        }
      }
      return num_added;
    }

    //! @brief sorts the population by fitness, then resizes it to the normal size
    template< typename t_MemberType, typename t_FitnessType>
    void EvolutionPopulation< t_MemberType, t_FitnessType>::Prune()
    {
      SortMembers();
      if( m_NormalPopulationSize > 0 && m_Members.GetSize() > m_NormalPopulationSize)
      {
        BCL_MessageVrb( "Resizing population to " + util::Format()( m_NormalPopulationSize) + " from " + util::Format()( m_Members.GetSize()));
        m_Members.Resize( m_NormalPopulationSize);
        RecalculateFitnessStatistics();
      }
    }

    //! @brief set the uniqueness criteria (used when AddMember is called)
    //! @param MEASURE the uniqueness measure to use
    template< typename t_MemberType, typename t_FitnessType>
    void EvolutionPopulation< t_MemberType, t_FitnessType>::SetUniquenessMeasure
    (
      const EvolutionMemberUniqueInterface< t_MemberType, t_FitnessType> &MEASURE
    )
    {
      m_MemberUnique = util::Implementation< EvolutionMemberUniqueInterface< t_MemberType, t_FitnessType> >( MEASURE);
    }

    //! @brief gets the uniqueness measure as an implementation
    //! @return the uniqueness measure
    template< typename t_MemberType, typename t_FitnessType>
    const util::Implementation< EvolutionMemberUniqueInterface< t_MemberType, t_FitnessType> >
      &EvolutionPopulation< t_MemberType, t_FitnessType>::GetUniquenessImplementation() const
    {
      return m_MemberUnique;
    }

    //! @brief get the uniqueness measure
    //! @return the uniqueness measure
    template< typename t_MemberType, typename t_FitnessType>
    const EvolutionMemberUniqueInterface< t_MemberType, t_FitnessType>
      &EvolutionPopulation< t_MemberType, t_FitnessType>::GetUniquenessMeasure() const
    {
      BCL_Assert
      (
        m_MemberUnique.IsDefined(),
        "GetUniquenessMeasure called, but uniqueness measure was unspecified for this population."
      );
      return *m_MemberUnique;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief assignment operator
    template< typename t_MemberType, typename t_FitnessType>
    EvolutionPopulation< t_MemberType, t_FitnessType> &EvolutionPopulation< t_MemberType, t_FitnessType>::operator =
    (
      const EvolutionPopulation< t_MemberType, t_FitnessType> &OTHER
    )
    {
      if( &OTHER == this)
      {
        return *this;
      }

      m_NormalPopulationSize = OTHER.m_NormalPopulationSize;
      m_MemberUnique = OTHER.m_MemberUnique;

      size_t final_size( OTHER.m_Members.GetSize());
      m_Members.AllocateMemory( final_size);

      for( size_t i( 0); i < final_size; ++i)
      {
        m_Members.PushBack
        ( 
          util::ShPtr< EvolutionPopulationMember< t_MemberType, t_FitnessType> >( OTHER.m_Members( i)->Clone())
        );
      }
      RegenerateFitnessVector();
      m_FitnessVectorStale = false;
      RecalculateFitnessStatistics();
      return *this;
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    template< typename t_MemberType, typename t_FitnessType>
    std::istream &EvolutionPopulation< t_MemberType, t_FitnessType>::Read( std::istream &ISTREAM)
    {
      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    template< typename t_MemberType, typename t_FitnessType>
    std::ostream &EvolutionPopulation< t_MemberType, t_FitnessType>::Write
    (
      std::ostream &OSTREAM, const size_t INDENT
    ) const
    {
      // end
      return OSTREAM;
    }

    //! @brief operator() that compares two members based on their fitnesses
    //! @param LHS the first member to compare
    //! @param RHS the second member to compare
    //! @return true if LHS is less than RHS
    template< typename t_MemberType, typename t_FitnessType>
    bool EvolutionMemberCompareLess< t_MemberType, t_FitnessType>::operator()
    (
      const EvolutionPopulationMember< t_MemberType, t_FitnessType> &LHS,
      const EvolutionPopulationMember< t_MemberType, t_FitnessType> &RHS
    ) const
    {
      return LHS.GetFitness() < RHS.GetFitness();
    }

    //! @brief operator() that compares two pointers to members based on their fitnesses
    //! @param LHS the first member to compare
    //! @param RHS the second member to compare
    //! @return true if value of LHS is less than value of RHS
    template< typename t_MemberType, typename t_FitnessType>
    bool EvolutionMemberCompareLess< t_MemberType, t_FitnessType>::operator()
    (
      const util::PtrInterface< EvolutionPopulationMember< t_MemberType, t_FitnessType> > &LHS,
      const util::PtrInterface< EvolutionPopulationMember< t_MemberType, t_FitnessType> > &RHS
    ) const
    {
      return operator()( *LHS, *RHS);
    }

  } // namespace opti
} // namespace bcl

#endif // BCL_OPTI_EVOLUTION_POPULATION_HPP_
