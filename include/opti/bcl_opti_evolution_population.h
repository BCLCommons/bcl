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

#ifndef BCL_OPTI_EVOLUTION_POPULATION_H_
#define BCL_OPTI_EVOLUTION_POPULATION_H_

// include the namespace header
#include "bcl_opti.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_opti_evolution_population_member.h"
#include "iterate/bcl_iterate_generic.h"
#include "linal/bcl_linal_vector.h"
#include "math/bcl_math_comparisons.h"
#include "math/bcl_math_discrete_set_selector.h"
#include "math/bcl_math_running_average_sd.h"
#include "util/bcl_util_implementation.h"
#include "util/bcl_util_serializable_interface.h"
#include "util/bcl_util_sh_ptr_vector.h"
#include "util/bcl_util_si_ptr_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opti
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class EvolutionPopulation
    //! @brief a class representing a population for use in an evolutionary algorithm.
    //! @details Has functions for inserting and selecting members, and calculates statistics on member fitnesses
    //!
    //! @see @link example_opti_evolution_population.cpp @endlink
    //! @author geanesar
    //! @date Oct 31, 2014
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_MemberType, typename t_FitnessType>
    class EvolutionPopulation :
      public util::ObjectInterface
    {
    private:

      //! The size of a "full" population
      size_t m_NormalPopulationSize;

      //! Members of the population
      util::ShPtrVector< EvolutionPopulationMember< t_MemberType, t_FitnessType> > m_Members;

      //! Statistical analyzer for fitnesses
      math::RunningAverageSD< t_FitnessType> m_FitnessAvg;

      //! vector containing fitnesses for rapid calculations (generated only when needed)
      mutable linal::Vector< t_FitnessType> m_FitnessVector;

      //! Whether the m_FitnessVector vector needs to be regenerated, e.g. after a compound was added
      mutable bool m_FitnessVectorStale;

      //! Pointer to a UniqueInterface for determining uniqueness
      util::Implementation< EvolutionMemberUniqueInterface< t_MemberType, t_FitnessType> > m_MemberUnique;

    public:

      //! @brief default constructor; creates population with unspecified normal size
      EvolutionPopulation();

      //! @brief constructor with population size
      EvolutionPopulation( const size_t &POPULATION_SIZE);

      //! @brief Constructor with population size and initial memberes
      //! @param POPULATION_SIZE the desired size of the population
      //! @param INITIAL_MEMBERS the members to add immediately
      EvolutionPopulation
      (
        const size_t &POPULATION_SIZE,
        util::ShPtrVector< EvolutionPopulationMember< t_MemberType, t_FitnessType> > &INITIAL_MEMBERS
      );

      //! @brief copy constructor, does a deep copy of the members
      EvolutionPopulation
      (
        const EvolutionPopulation< t_MemberType, t_FitnessType> &OTHER
      );

      //! @brief copies this class
      //! @return A copy of this class
      EvolutionPopulation *Clone() const;

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief Get the members from the population
      //! @return Returns a const reference to all members
      util::SiPtrVector< const EvolutionPopulationMember< t_MemberType, t_FitnessType> > GetMembers() const;

      //! @brief gets an iterator for members
      //! @return a generic iterator to members of the populations
      iterate::Generic< const EvolutionPopulationMember< t_MemberType, t_FitnessType> > GetMembersIterator() const;

      //! @brief gets an iterator for members
      //! @return a generic iterator to members of the populations
      iterate::Generic< EvolutionPopulationMember< t_MemberType, t_FitnessType> > GetMembersIterator();

      //! @brief Return the top percentage of fittest members
      //! @param PERCENTAGE the percentage of members to return, starting with the fittest
      //! @return SiPtrVector containing pairs of members and their fitness scores
      util::SiPtrVector< const EvolutionPopulationMember< t_MemberType, t_FitnessType> > GetFittestPercent
      (
        const float &PERCENTAGE
      ) const;

      //! @brief Get a number of the fittest members of a population
      //! @param NUMBER the number of members to return, starting with the fittest
      //! @return SiPtrVector to fittest members
      util::SiPtrVector< const EvolutionPopulationMember< t_MemberType, t_FitnessType> > GetFittestMembers
      (
        const size_t &NUMBER
      ) const;

      //! @brief get members with a fitness over a cutoff
      //! @param FITNESS minimum fitness
      //! @return SiPtrVector
      util::SiPtrVector< const EvolutionPopulationMember< t_MemberType, t_FitnessType> > GetMembersWithMinimumFitness
      (
        const t_FitnessType &FITNESS
      ) const;

      //! @brief return a member using a probability distribution
      //! @details this copies the set selector, so if this is expensive then there may be a better way to do this
      //! @param SELECTOR the discrete set selector (i.e. probability distribution function) to use
      //! @return SiPtr to the chosen member
      util::SiPtr< const EvolutionPopulationMember< t_MemberType, t_FitnessType> > GetRandomMemberFromDistribution
      (
        const math::DiscreteSetSelector &SELECTOR
      ) const;

      //! @brief Get a random member of the population (uniform probability)
      //! @return a pointer to the chosen member
      util::SiPtr< const EvolutionPopulationMember< t_MemberType, t_FitnessType> > GetRandomMember() const;

      //! @brief Get the current size of the population
      //! @return The number of members currently held in the population
      size_t GetCurrentSize() const;

      //! @brief Get the normal size of the population, i.e. the size of a full population
      //! @return the size of a full population
      size_t GetNormalSize() const;

      //! @brief set the normal size of the population
      //! @param SIZE the size
      void SetNormalSize( const size_t &SIZE);

      //! @brief Clears the population (deletes members)
      void Reset();

      //! @brief Check if the population has a particular member, then add it
      //! @param MEMBER the member to add
      //! @return true if the member was added, false otherwise
      bool HasMember( const EvolutionPopulationMember< t_MemberType, t_FitnessType> &MEMBER) const;

      //! @brief Gets the lowest fitness present in the population
      //! @return The value of the lowest fitness in the population, or undefined if the population is empty
      t_FitnessType GetLowestFitness() const;

      //! @brief Gets the highest fitness in the population
      //! @return The value of the highest fitness in the population, or undefined if the population is empty
      t_FitnessType GetHighestFitness() const;

      //! @brief Gets the average fitness in the population
      //! @return The value of the average fitness of the population
      t_FitnessType GetFitnessAverage() const;

      //! @brief Gets the standard deviation of the fitnesses
      //! @return the fitness standard deviation
      t_FitnessType GetFitnessStdDev() const;

      //! @brief Adds a member to the population
      //! @param MEMBER the member to add
      //! @param FITNESS_SCORE the fitness of the new member
      //! @return true if the member was added successfully
      bool AddMember( const t_MemberType &MEMBER, const t_FitnessType &FITNESS_SCORE);

      //! @brief add a member to the population, then re-sort the members
      //! @param NEW_MEMBER the member to add
      bool AddMember( const EvolutionPopulationMember< t_MemberType, t_FitnessType> &NEW_MEMBER);

      //! @brief add a member to the population
      //! @param NEW_MEMBER the member to add
      bool AddMember( const util::ShPtr< EvolutionPopulationMember< t_MemberType, t_FitnessType> > &NEW_MEMBER);

      //! @brief add members to the population
      //! @param MEMBERS the members to add
      size_t AddMembers( const util::ShPtrVector< EvolutionPopulationMember< t_MemberType, t_FitnessType> > &MEMBERS);

      //! @brief sorts the population by fitness, then resizes it to the normal size
      void Prune();

      //! @brief set the uniqueness criteria (used when AddMember is called)
      //! @param MEASURE the uniqueness measure to use
      void SetUniquenessMeasure( const EvolutionMemberUniqueInterface< t_MemberType, t_FitnessType> &MEASURE);

      //! @brief gets the uniqueness measure as an implementation
      //! @return the uniqueness measure
      const util::Implementation< EvolutionMemberUniqueInterface< t_MemberType, t_FitnessType> >
        &GetUniquenessImplementation() const;

      //! @brief get the uniqueness measure
      //! @return the uniqueness measure
      const EvolutionMemberUniqueInterface< t_MemberType, t_FitnessType> &GetUniquenessMeasure() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief assignment operator
      EvolutionPopulation< t_MemberType, t_FitnessType> &operator =( const EvolutionPopulation< t_MemberType, t_FitnessType> &OTHER);

    private:

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief sorts members by fitness score from highest to lowest
      void SortMembers();

      //! @brief Recalculates the fitness statistics (average, std dev)
      void RecalculateFitnessStatistics();

      //! @brief regenerates m_FitnessVector from current members
      void RegenerateFitnessVector() const;

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    };

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class EvolutionMemberUniqueInterface
    //! @brief virtual comparison class for determining if an population member is already included in a population
    //!
    //! @remarks example unnecessary
    //! @author geanesar
    //! @date Oct 31, 2014
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_MemberType, typename t_FitnessType>
    class EvolutionMemberUniqueInterface :
      public util::SerializableInterface
    {
    public:

      //! @brief Clone function
      //! @return a pointer to a new copy of the class (must be overridden)
      virtual EvolutionMemberUniqueInterface *Clone() const = 0;

      //! @brief function for testing if a member is contained in the population
      //! @param MEMBER the member to check
      //! @param POPULATIon the population to search
      //! @return true if the member is contained in the population, false otherwise
      virtual bool Contains
      (
        const EvolutionPopulationMember< t_MemberType, t_FitnessType> &MEMBER,
        const EvolutionPopulation< t_MemberType, t_FitnessType> &POPULATION
      ) const = 0;

      //! @brief operator for comparing members.  An alias for Contains()
      //! @param MEMBER the member to look for
      //! @param POPULATION the population to search in
      //! @return true if the population contains the member, false otherwise
      bool operator()
      (
        const EvolutionPopulationMember< t_MemberType, t_FitnessType> &MEMBER,
        const EvolutionPopulation< t_MemberType, t_FitnessType> &POPULATION
      ) const
      {
        return Contains( MEMBER, POPULATION);
      }

    };

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class EvolutionMemberCompareLess
    //! @brief A binary comparison class EvolutionPopulationMembers; compares their fitness scores
    //!
    //! @remarks example unnecessary
    //! @author geanesar
    //! @date Oct 31, 2014
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_MemberType, typename t_FitnessType>
    struct EvolutionMemberCompareLess
    {
      //! @brief operator() that compares two members based on their fitnesses
      //! @param LHS the first member to compare
      //! @param RHS the second member to compare
      //! @return if LHS is less than RHS
      bool operator()
      (
        const EvolutionPopulationMember< t_MemberType, t_FitnessType> &LHS,
        const EvolutionPopulationMember< t_MemberType, t_FitnessType> &RHS
      ) const;

      //! @brief operator() that compares two pointers to members based on their fitnesses
      //! @param LHS the first member to compare
      //! @param RHS the second member to compare
      //! @return if value of LHS is less than value of RHS
      bool operator()
      (
        const util::PtrInterface< EvolutionPopulationMember< t_MemberType, t_FitnessType> > &LHS,
        const util::PtrInterface< EvolutionPopulationMember< t_MemberType, t_FitnessType> > &RHS
      ) const;

    };

  } // namespace opti
} // namespace bcl

// include template implementation file
#include "bcl_opti_evolution_population.hpp"

#endif // BCL_OPTI_EVOLUTION_POPULATION_H_
