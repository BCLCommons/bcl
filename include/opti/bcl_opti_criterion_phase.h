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

#ifndef BCL_OPTI_CRITERION_PHASE_H_
#define BCL_OPTI_CRITERION_PHASE_H_

// include the namespace header
#include "bcl_opti.h"

// includes from bcl - sorted alphabetically
#include "bcl_opti_criterion_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opti
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class CriterionPhase
    //! @brief This criteria is met if the tracker has the specified phase.
    //!
    //! @tparam t_ArgumentType is the type of the approximation argument
    //! @tparam t_ResultType is the type of the approximation result
    //!
    //! @see @link example_opti_criterion_phase.cpp @endlink
    //! @author fischea
    //! @date Dec 14, 2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_ArgumentType, typename t_ResultType>
    class CriterionPhase :
      public CriterionInterface< t_ArgumentType, t_ResultType>
    {

    //////////
    // data //
    //////////

    private:

      //! the phases for which the criterion is met
      storage::Set< PhaseEnum> m_Phases;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // @brief default constructor
      CriterionPhase() :
        m_Phases()
      {
      }

      //! @brief construct from maximum number of allowed iterations
      //! @param PHASES the phases for which the criterion is met
      CriterionPhase( const storage::Set< PhaseEnum> &PHASES) :
        m_Phases( PHASES)
      {
      }

      //! @brief Clone function
      //! @return pointer to a new CriterionPhase< t_ArgumentType, t_ResultType>
      CriterionPhase *Clone() const
      {
        return new CriterionPhase< t_ArgumentType, t_ResultType>( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns the class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief returns whether the maximum number of iterations has been observed for the given tracker
      //! @param TRACKER tracker to evaluate the criterion for
      //! @return true, if any of the termination criteria are met yet
      bool CriteriaMet( const Tracker< t_ArgumentType, t_ResultType> &TRACKER) const
      {
        // check if the phase of the tracker is equals to the criterion phases
        for
        (
          storage::Set< PhaseEnum>::const_iterator it( m_Phases.Begin()), it_end( m_Phases.End());
          it != it_end;
          ++it
        )
        {
          if( PhasesEqual( *it, TRACKER.GetPhase()))
          {
            return true;
          }
        }

        return false;
      }

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const
      {
        io::Serializer serializer;
        serializer.SetClassDescription( "Triggers at a specified phase");
        serializer.AddInitializer
        (
          "",
          "",
          io::Serialization::GetAgent( &m_Phases)
        );
        return serializer;
      }

    }; // template class CriterionPhase< t_ArgumentType, t_ResultType>

    // instantiate s_Instance
    template< typename t_ArgumentType, typename t_ResultType>
    const util::SiPtr< const util::ObjectInterface> CriterionPhase< t_ArgumentType, t_ResultType>::s_Instance
    (
      util::Enumerated< CriterionInterface< t_ArgumentType, t_ResultType> >::AddInstance( new CriterionPhase< t_ArgumentType, t_ResultType>())
    );

  } // namespace opti
} // namespace bcl

#endif // BCL_OPTI_CRITERION_PHASE_H_
