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

#ifndef BCL_OPTI_CRITERION_IMPROVEMENT_RATIO_H_
#define BCL_OPTI_CRITERION_IMPROVEMENT_RATIO_H_

// include the namespace header
#include "bcl_opti.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_opti_criterion_interface.h"
#include "io/bcl_io_serialization.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opti
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class CriterionImprovementRatio
    //! @brief CriterionInterface implementation that uses a decider function to evaluate whether the criterion is met
    //!
    //! @tparam t_ArgumentType is the type of the approximation argument
    //! @tparam t_ResultType is the type of the approximation result
    //!
    //! @see @link example_opti_criterion_improvement_ratio.cpp @endlink
    //! @author karakam, fischea, mendenjl
    //! @date Sep 03, 2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_ArgumentType, typename t_ResultType>
    class CriterionImprovementRatio :
      public CriterionInterface< t_ArgumentType, t_ResultType>
    {

    //////////
    // data //
    //////////

    private:

      //! improvement ratio required
      double m_Ratio;

      //! number of steps within this improvement ratio has to occur
      size_t m_NumberSteps;

      //! last iteration with improvement above ratio
      mutable size_t m_LastSatisfactoryIteration;

      //! Last satisfactory result
      mutable t_ResultType m_LastSatisfactoryResult;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // @brief default constructor
      CriterionImprovementRatio() :
        m_Ratio( 0.0),
        m_NumberSteps( 0),
        m_LastSatisfactoryIteration( util::GetUndefined< size_t>()),
        m_LastSatisfactoryResult()
      {
      }

      //! @brief construct from required improvement ratio and step number
      //! @param RATIO required improvement ratio
      //! @param NUMBER_STEPS number of steps within the improvement ratio has to occur
      CriterionImprovementRatio( const double &RATIO, const size_t &NUMBER_STEPS) :
        m_Ratio( 0.0),
        m_NumberSteps( NUMBER_STEPS),
        m_LastSatisfactoryIteration( util::GetUndefined< size_t>()),
        m_LastSatisfactoryResult()
      {
      }

      //! @brief Clone function
      //! @return pointer to a new CriterionImprovementRatio< t_ArgumentType, t_ResultType>
      CriterionImprovementRatio *Clone() const
      {
        return new CriterionImprovementRatio< t_ArgumentType, t_ResultType>( *this);
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
        static const std::string s_alias( "ImprovementRatio");
        return s_alias;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief returns the function value of the decider function the the current model
      //! @param TRACKER tracker to evaluate the criterion for
      //! @return true, if any of the termination criteria are met yet
      bool CriteriaMet( const Tracker< t_ArgumentType, t_ResultType> &TRACKER) const
      {
        // if no result has been accepted yet
        if( !util::IsDefined( m_LastSatisfactoryIteration))
        {
          m_LastSatisfactoryIteration = TRACKER.GetIteration();
          m_LastSatisfactoryResult = TRACKER.GetCurrent()->Second();

          return false;
        }

        // calculate the ratio of improvement
        const double ratio
        (
          ( m_LastSatisfactoryResult - TRACKER.GetCurrent()->Second()) / m_LastSatisfactoryResult
        );

        if( ratio >= m_Ratio && TRACKER.GetIteration() - m_LastSatisfactoryIteration < m_NumberSteps)
        {
          m_LastSatisfactoryIteration = TRACKER.GetIteration();
          m_LastSatisfactoryResult = TRACKER.GetCurrent()->Second();

          return false;
        }

        // return true if there are still more steps possible
        return TRACKER.GetIteration() - m_LastSatisfactoryIteration >= m_NumberSteps;
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
        serializer.SetClassDescription
        (
          "Criteria that the ratio of improved results is at least X within Y iterations"
        );
        serializer.AddInitializer
        (
          "max iterations",
          "Max iterations allowed to reach the desired ratio of result",
          io::Serialization::GetAgent( &m_NumberSteps)
        );
        serializer.AddInitializer
        (
          "ratio",
          "Ratio of results must be at least this much within the given iterations",
          io::Serialization::GetAgent( &m_Ratio)
        );
        serializer.AddDataMember
        (
          "last satisfactory iteration",
          io::Serialization::GetAgent( &m_LastSatisfactoryIteration)
        );
        serializer.AddDataMember
        (
          "last satisfactory result",
          io::Serialization::GetAgent( &m_LastSatisfactoryResult)
        );
        return serializer;
      }

    }; // template class CriterionImprovementRatio< t_ArgumentType, t_ResultType>

    // instantiate s_Instance
    template< typename t_ArgumentType, typename t_ResultType>
    const util::SiPtr< const util::ObjectInterface> CriterionImprovementRatio< t_ArgumentType, t_ResultType>::s_Instance
    (
      util::Enumerated< CriterionInterface< t_ArgumentType, t_ResultType> >::AddInstance( new CriterionImprovementRatio< t_ArgumentType, t_ResultType>())
    );

  } // namespace opti
} // namespace bcl

#endif // BCL_OPTI_CRITERION_IMPROVEMENT_RATIO_H_
