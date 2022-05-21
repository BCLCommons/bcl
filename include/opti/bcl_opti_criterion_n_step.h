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

#ifndef BCL_OPTI_CRITERION_N_STEP_H_
#define BCL_OPTI_CRITERION_N_STEP_H_

// include the namespace header
#include "bcl_opti.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_opti_criterion_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opti
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class CriterionNStep
    //! @brief CriterionInterface implementation that is met very n-th step.
    //!
    //! @tparam t_ArgumentType is the type of the approximation argument
    //! @tparam t_ResultType is the type of the approximation result
    //!
    //! @see @link example_opti_criterion_n_step.cpp @endlink
    //! @author fischea
    //! @date Jan 25, 2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_ArgumentType, typename t_ResultType>
    class CriterionNStep :
      public CriterionInterface< t_ArgumentType, t_ResultType>
    {

    //////////
    // data //
    //////////

    private:

      //! repeat interval
      mutable size_t m_RepeatInterval;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      CriterionNStep() :
        m_RepeatInterval( 0)
      {
      }

      //! @brief construct from repeat interval
      //! @param REPEAT_INTERVAL maximum number of repeats
      CriterionNStep( const size_t &REPEAT_INTERVAL) :
        m_RepeatInterval( REPEAT_INTERVAL)
      {
      }

      //! @brief Clone function
      //! @return pointer to a new CriterionNStep< t_ArgumentType, t_ResultType>
      CriterionNStep *Clone() const
      {
        return new CriterionNStep< t_ArgumentType, t_ResultType>( *this);
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

      //! @brief returns the class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetAlias() const
      {
        static const std::string s_alias( "EveryNSteps");
        return s_alias;
      }

      //! @brief returns the maximum number of repeats within tolerance
      //! @return maximum number of repeats within tolerance
      size_t GetRepeatInterval() const
      {
        return m_RepeatInterval;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief returns whether the termination criterion is met
      //! @param TRACKER tracker to evaluate the criterion for
      //! @return true, if the criterion is met
      bool CriteriaMet( const Tracker< t_ArgumentType, t_ResultType> &TRACKER) const
      {
        return TRACKER.GetIteration() % m_RepeatInterval == 0;
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
        serializer.SetClassDescription( "Triggers after every N iterations");
        serializer.AddInitializer
        (
          "",
          "",
          io::Serialization::GetAgent( &m_RepeatInterval)
        );
        return serializer;
      }

    }; // template class CriterionNStep< t_ArgumentType, t_ResultType>

    // instantiate s_Instance
    template< typename t_ArgumentType, typename t_ResultType>
    const util::SiPtr< const util::ObjectInterface> CriterionNStep< t_ArgumentType, t_ResultType>::s_Instance
    (
      util::Enumerated< CriterionInterface< t_ArgumentType, t_ResultType> >::AddInstance( new CriterionNStep< t_ArgumentType, t_ResultType>())
    );

  } // namespace opti
} // namespace bcl

#endif // BCL_OPTI_CRITERION_N_STEP_H_
