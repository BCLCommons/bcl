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

#ifndef BCL_OPTI_CRITERION_ELAPSED_TIME_H_
#define BCL_OPTI_CRITERION_ELAPSED_TIME_H_

// include the namespace header
#include "bcl_opti.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_opti_criterion_interface.h"
#include "command/bcl_command_parameter_check_ranged.h"
#include "io/bcl_io_serialization_via_static_functions.h"
#include "util/bcl_util_stopwatch.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opti
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class CriterionElapsedTime
    //! @brief Termination criterion for approximator framework that triggers after specified amount of time.
    //! @details This criterion is useful for letting models train a defined duration. problems with job cancellation
    //! on a cluster due to overtime can be avoided this way.
    //!
    //! @tparam t_ArgumentType is the type of the approximation argument
    //! @tparam t_ResultType is the type of the approximation result
    //!
    //! @see @link example_opti_criterion_elapsed_time.cpp @endlink
    //! @author butkiem1, fischea
    //! @date Jan 4, 2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_ArgumentType, typename t_ResultType>
    class CriterionElapsedTime :
      public CriterionInterface< t_ArgumentType, t_ResultType>
    {

    //////////
    // data //
    //////////

    private:

      //! duration time until the terminate triggers
      util::Time m_Duration;

      //! stop watch that measures duration time
      mutable util::Stopwatch m_Timer;

    public:

      //! stop watch descriptor
      static const std::string s_StopWatchDescriptor;

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      CriterionElapsedTime() :
        m_Duration(),
        m_Timer( s_StopWatchDescriptor)
      {
        m_Timer.Start();
      }

      //! @brief construct from maximum duration time
      //! @param DURATION_TIME maximum duration time
      CriterionElapsedTime( const util::Time &DURATION_TIME) :
        m_Duration( DURATION_TIME),
        m_Timer( util::Stopwatch( s_StopWatchDescriptor))
      {
        m_Timer.Start();
      }

      //! @brief Clone function
      //! @return pointer to a new CriterionElapsedTime< t_ArgumentType, t_ResultType>
      CriterionElapsedTime *Clone() const
      {
        return new CriterionElapsedTime< t_ArgumentType, t_ResultType>( *this);
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
        static const std::string s_alias( "Minutes");
        return s_alias;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief returns if any of the combined termination criteria are met
      //! @param RACKER tracker to evaluate the criterion for
      //! @return true, if any of the termination criteria are met
      bool CriteriaMet( const Tracker< t_ArgumentType, t_ResultType> &TRACKER) const
      {
        // check whether total stop watch time exceeds duration
        return ( m_Timer.GetTotalTime() > m_Duration);
      }

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief hook that derived classes can override, called after TryRead successfully reads a serializer containing only initializer info (no data variables)
      //! @param SERIALIZER the serializer object with initialization information
      //! @param ERR_STREAM stream to write out errors to
      //! @return true, unless there were new errors
      bool ReadInitializerSuccessHook( const util::ObjectDataLabel &SERIALIZER, std::ostream &ERR_STREAM)
      {
        m_Timer.Reset();
        m_Timer.Start();
        return true;
      }

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const
      {
        io::Serializer serializer;
        serializer.SetClassDescription
        (
          "Maximum # of minutes for optimization"
        );
        serializer.AddInitializer
        (
          "",
          "",
          util::OwnPtr< io::SerializationBase< util::Time> >
          (
            new io::SerializationViaStaticFunctions< util::Time>
            (
              command::ParameterCheckRanged< double>( 0.0, 60.0 * 24.0),
              &util::Time::ConvertTimeToMinutesString,
              &util::Time::CreateTimeFromMinutesString,
              &m_Duration
            )
          )
        );
        return serializer;
      }

    }; // template class CriterionElapsedTime< t_ArgumentType, t_ResultType>

    // instantiate s_Instance
    template< typename t_ArgumentType, typename t_ResultType>
    const util::SiPtr< const util::ObjectInterface> CriterionElapsedTime< t_ArgumentType, t_ResultType>::s_Instance
    (
      util::Enumerated< CriterionInterface< t_ArgumentType, t_ResultType> >::AddInstance( new CriterionElapsedTime< t_ArgumentType, t_ResultType>())
    );

    // initialize descriptor for stop watch
    template< typename t_ArgumentType, typename t_ResultType>
    const std::string CriterionElapsedTime< t_ArgumentType, t_ResultType>::s_StopWatchDescriptor( "Timer for TerminateOnElapsedTime");

  } // namespace opti
} // namespace bcl

#endif // BCL_OPTI_CRITERION_ELAPSED_TIME_H_
