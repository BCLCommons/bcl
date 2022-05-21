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

#ifndef BCL_OPTI_EVOLUTION_POPULATION_MEMBER_H_
#define BCL_OPTI_EVOLUTION_POPULATION_MEMBER_H_

// include the namespace header
#include "bcl_opti.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialize.h"
#include "math/bcl_math_binary_function_interface.h"
#include "util/bcl_util_binary_function_interface.h"
#include "util/bcl_util_object_data_label.h"
#include "util/bcl_util_object_interface.h"

namespace bcl
{
  namespace opti
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class EvolutionPopulationMember
    //! @brief a container class for members of an EvolutionPopulation
    //! @details This class is a container for information relevant to members of an EvolutionPopulation
    //!
    //! @see @link example_opti_evolution_population_member.cpp @endlink
    //! @author geanesar
    //! @date Oct 27 2014
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_MemberType, typename t_FitnessType>
    class EvolutionPopulationMember :
      public util::ObjectInterface
    {
    private:

      //! the member
      t_MemberType m_Member;

      //! the member's fitness
      t_FitnessType m_Fitness;

      //! the member's history
      util::ObjectDataLabel m_MemberHistory;

    public:

      //! @brief default constructor
      EvolutionPopulationMember() :
        m_Member(),
        m_Fitness( util::GetUndefined< t_FitnessType>()),
        m_MemberHistory()
      {
      }

      //! @breif constructor
      //! @param MEMBER the member to store
      //! @param FITNESS the member's fitness
      EvolutionPopulationMember
      (
        const t_MemberType &MEMBER,
        const t_FitnessType &FITNESS
      ) :
        m_Member( MEMBER),
        m_Fitness( FITNESS)
      {
      }

      //! @breif constructor
      //! @param MEMBER the member to store
      //! @param FITNESS the member's fitness
      //! @param HISTORY the member's history
      EvolutionPopulationMember
      (
        const t_MemberType &MEMBER,
        const t_FitnessType &FITNESS,
        const util::ObjectDataLabel &HISTORY
      ) :
        m_Member( MEMBER),
        m_Fitness( FITNESS),
        m_MemberHistory( HISTORY)
      {
      }

      //! @brief clone operation
      //! @return a pointer to a new copy of this class
      EvolutionPopulationMember *Clone() const
      {
        return new EvolutionPopulationMember( *this);
      }

      //! @brief returns class name of this object
      //! @return the class name
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief get the name of this object when used in a dynamic context
      //! @return the name of this object when used in a dynamic context
      const t_MemberType &GetMember() const
      {
        return m_Member;
      }

      //! @brief set the member data
      //! @param MEMBER the member
      void SetMember( const t_MemberType &MEMBER)
      {
        m_Member = MEMBER;
      }

      //! @brief Get the fitness of the member
      //! @return the member's fitness
      const t_FitnessType &GetFitness() const
      {
        return m_Fitness;
      }

      //! @brief Set fitness score
      //! @param FITNESS the fitness score
      void SetFitness( const t_FitnessType &FITNESS)
      {
        m_Fitness = FITNESS;
      }

      //! @brief Get the history of the member
      //! @return the member's history
      const util::ObjectDataLabel &GetHistory() const
      {
        return m_MemberHistory;
      }

      //! @brief Set the history for the member
      //! @param HISTORY the history
      void SetHistory( const util::ObjectDataLabel &HISTORY)
      {
        m_MemberHistory = HISTORY;
      }

    ///////////////
    // operators //
    ///////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM)
      {
        io::Serialize::Read( m_Member, ISTREAM);
        io::Serialize::Read( m_Fitness, ISTREAM);
        io::Serialize::Read( m_MemberHistory, ISTREAM);
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT - number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        io::Serialize::Write( m_Member, OSTREAM, INDENT);
        io::Serialize::Write( m_Fitness, OSTREAM, INDENT);
        io::Serialize::Write( m_MemberHistory, OSTREAM, INDENT);
        return OSTREAM;
      }
    };

  } // namespace opti
} // namespace bcl

#endif // BCL_OPTI_EVOLUTION_POPULATION_MEMBER_H_
