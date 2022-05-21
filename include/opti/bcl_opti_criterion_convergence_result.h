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

#ifndef BCL_OPTI_CRITERION_CONVERGENCE_RESULT_H_
#define BCL_OPTI_CRITERION_CONVERGENCE_RESULT_H_

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
    //! @class CriterionConvergenceResult
    //! @brief CriterionInterface implementation for terminating after a repetitive difference between a previous and
    //! current result within a given tolerance
    //!
    //! @tparam t_ArgumentType is the type of the approximation argument.
    //! @tparam t_ResultType is the type of the approximation result. A difference has to be supported and a '<'
    //! less than operation
    //!
    //! @see @link example_opti_criterion_convergence_result.cpp @endlink
    //! @author woetzen, fischea
    //! @date Dec 13, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_ArgumentType, typename t_ResultType>
    class CriterionConvergenceResult :
      public CriterionInterface< t_ArgumentType, t_ResultType>
    {

    //////////
    // data //
    //////////

    private:

      //! counts the number of repeats
      mutable size_t m_NumberRepeats;

      //! maximum number of repeats
      size_t m_MaxNumberRepeats;

      //! tolerance under which arguments are considered as repeats
      t_ResultType m_Tolerance;

      //! previous result
      mutable t_ResultType m_PreviousResult;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

      //! @brief default constructor
      CriterionConvergenceResult() :
        m_NumberRepeats( util::GetUndefined< size_t>()),
        m_MaxNumberRepeats( 0),
        m_Tolerance( 0.0)
      {
      }

      //! @brief construct from maximum number of repeats and tolerance
      //! @param MAX_NUMBER_REPEATS maximum number of repeats
      //! @param TOLERANCE tolerance to consider results as repeats
      CriterionConvergenceResult
      (
        const size_t &MAX_NUMBER_REPEATS,
        const t_ResultType &TOLERANCE
      ) :
        m_NumberRepeats( util::GetUndefined< size_t>()),
        m_MaxNumberRepeats( MAX_NUMBER_REPEATS),
        m_Tolerance( TOLERANCE)
      {
      }

      //! @brief Clone function
      //! @return pointer to a new CriterionConvergenceResult< t_ArgumentType, t_ResultType>
      CriterionConvergenceResult *Clone() const
      {
        return new CriterionConvergenceResult< t_ArgumentType, t_ResultType>( *this);
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
        static const std::string s_alias( "ConvergenceResult");
        return s_alias;
      }

      //! @brief returns the maximum number of repeats within tolerance
      //! @return maximum number of repeats within tolerance
      size_t GetMaxNumberRepeats() const
      {
        return m_MaxNumberRepeats;
      }

      //! @brief returns the current number of repeats within tolerance
      //! @return current number of repeats within tolerance
      size_t GetNumberRepeats() const
      {
        return m_NumberRepeats;
      }

      //! @brief returns the tolerance to consider arguments as repeats
      //! @return tolerance to consider arguments as repeats
      t_ResultType GetTolerance() const
      {
        return m_Tolerance;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief returns whether the maximum number of iterations with difference of two successive results within the
      //! tolerance have been observed
      //! @param TRACKER tracker to evaluate the criterion for
      //! @return true, if any of the termination criteria are met yet
      bool CriteriaMet( const Tracker< t_ArgumentType, t_ResultType> &TRACKER) const
      {
        if( !TRACKER.GetCurrent().IsDefined())
        {
          return false;
        }

        // get the current result pair from the tracker
        const t_ResultType &current( TRACKER.GetCurrent()->Second());

        // if no previous argument result pair has been observed, the criterion is not met
        if( !util::IsDefined( m_NumberRepeats))
        {
          m_PreviousResult = current;
          m_NumberRepeats = 0;
          return false;
        }

        // check whether the new argument lies within the tolerance
        if( math::EqualWithinTolerance( current, m_PreviousResult, m_Tolerance))
        {
          ++m_NumberRepeats;
        }
        else
        {
          m_NumberRepeats = 0;
        }

        // set previous tracker result to current item
        m_PreviousResult = current;

        // check whether criteria is met
        return m_NumberRepeats >= m_MaxNumberRepeats;
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
          "Criteria that the result converges to a steady value for a given number of iterations"
        );
        serializer.AddInitializer
        (
          "repeats",
          "Required # of consecutive iterations in which the result has not changed more than the tolerance",
          io::Serialization::GetAgent( &m_MaxNumberRepeats)
        );
        serializer.AddInitializer
        (
          "tolerance",
          "tolerance under which results are considered as repeats",
          io::Serialization::GetAgent( &m_Tolerance)
        );
        serializer.AddDataMember( "current repetitions", io::Serialization::GetAgent( &m_NumberRepeats));
        serializer.AddDataMember( "last result", io::Serialization::GetAgent( &m_PreviousResult));
        return serializer;
      }

    }; // template class CriterionConvergenceResult< t_ArgumentType, t_ResultType>

    // instantiate s_Instance
    template< typename t_ArgumentType, typename t_ResultType>
    const util::SiPtr< const util::ObjectInterface> CriterionConvergenceResult< t_ArgumentType, t_ResultType>::s_Instance
    (
      util::Enumerated< CriterionInterface< t_ArgumentType, t_ResultType> >::AddInstance( new CriterionConvergenceResult< t_ArgumentType, t_ResultType>())
    );

  } // namespace opti
} // namespace bcl

#endif // BCL_OPTI_CRITERION_CONVERGENCE_RESULT_H_
