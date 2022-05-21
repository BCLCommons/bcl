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

#ifndef BCL_OPTI_CRITERION_RESULT_THRESHOLD_H_
#define BCL_OPTI_CRITERION_RESULT_THRESHOLD_H_

// include the namespace header
#include "bcl_opti.h"

// includes from bcl - sorted alphabetically
#include "bcl_opti_criterion_interface.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opti
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class CriterionResultThreshold
    //! @brief CriterionInterface implementation for terminating when the result is better than a particular threshold
    //!
    //! @tparam t_ArgumentType is the type of the approximation argument
    //! @tparam t_ResultType is the type of the approximation result
    //!
    //! @see @link example_opti_criterion_result_threshold.cpp @endlink
    //! @author mendenjl
    //! @date Aug 30, 2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_ArgumentType, typename t_ResultType>
    class CriterionResultThreshold :
      public CriterionInterface< t_ArgumentType, t_ResultType>
    {

    //////////
    // data //
    //////////

    private:

      //! Desired threshold value; once the result is better than this, the criterion triggers
      t_ResultType m_Threshold;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      CriterionResultThreshold() :
        m_Threshold()
      {
      }

      //! @brief construct from threshold
      //! @param THRESHOLD Thresold value; if the result improves on this, it triggers
      CriterionResultThreshold( const t_ResultType &THRESHOLD) :
        m_Threshold( THRESHOLD)
      {
      }

      //! @brief Clone function
      //! @return pointer to a new CriterionResultThreshold< t_ArgumentType, t_ResultType>
      CriterionResultThreshold *Clone() const
      {
        return new CriterionResultThreshold< t_ArgumentType, t_ResultType>( *this);
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
        static const std::string s_alias( "Threshold");
        return s_alias;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief returns if the termination criterion is met for the given tracker
      //! @param TRACKER tracker to evaluate the criterion for
      //! @return true, if any of the termination criteria are met yet
      bool CriteriaMet( const Tracker< t_ArgumentType, t_ResultType> &TRACKER) const
      {
        if( !TRACKER.GetBest().IsDefined())
        {
          return false;
        }
        return DoesImprove( TRACKER.GetBest()->Second(), m_Threshold, TRACKER.GetImprovementType());
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
        serializer.SetClassDescription( "Triggers when the result is better than the given value");
        serializer.AddInitializer
        (
          "",
          "",
          io::Serialization::GetAgent( &m_Threshold)
        );
        return serializer;
      }

    }; // template class CriterionResultThreshold< t_ArgumentType, t_ResultType>

    // instantiate s_Instance
    template< typename t_ArgumentType, typename t_ResultType>
    const util::SiPtr< const util::ObjectInterface> CriterionResultThreshold< t_ArgumentType, t_ResultType>::s_Instance
    (
      util::Enumerated< CriterionInterface< t_ArgumentType, t_ResultType> >::AddInstance( new CriterionResultThreshold< t_ArgumentType, t_ResultType>())
    );

  } // namespace opti
} // namespace bcl

#endif // BCL_OPTI_CRITERION_RESULT_THRESHOLD_H_
