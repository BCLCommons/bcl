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

#ifndef BCL_OPTI_CRITERION_RESULT_CHANGED_H_
#define BCL_OPTI_CRITERION_RESULT_CHANGED_H_

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
    //! @class CriterionResultChanged
    //! @brief CriterionInterface implementation for termination after a certain number of unimproved steps in a row.
    //!
    //! @tparam t_ArgumentType is the type of the approximation argument.
    //! @tparam t_ResultType is the type of the approximation result. A difference has to be supported and a '<' less
    //! than operation
    //!
    //! @see @link example_opti_criterion_result_changed.cpp @endlink
    //! @author mendenjl
    //! @date Jul 18, 2017
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_ArgumentType, typename t_ResultType>
    class CriterionResultChanged :
      public CriterionInterface< t_ArgumentType, t_ResultType>
    {

    //////////
    // data //
    //////////

      mutable t_ResultType m_LastResultValue;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

      //! @brief constructor
      CriterionResultChanged() :
        m_LastResultValue( util::GetUndefined< t_ResultType>())
      {
      }

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief Clone function
      //! @return pointer to a new CriterionResultChanged< t_ArgumentType, t_ResultType>
      CriterionResultChanged *Clone() const
      {
        return new CriterionResultChanged< t_ArgumentType, t_ResultType>( *this);
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
        static const std::string s_alias( "Unchanged");
        return s_alias;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief returns whether the maximum number of iterations with difference of two successive results within the
      //! tolerance has been observed
      //! @param TRACKER tracker to evaluate the criterion for
      //! @return true, if the maximum number of iterations with difference of two successive results within the
      //! tolerance has been observed
      bool CriteriaMet( const Tracker< t_ArgumentType, t_ResultType> &TRACKER) const
      {
        if( TRACKER.GetCurrent()->Second() == m_LastResultValue)
        {
          return false;
        }
        m_LastResultValue = TRACKER.GetCurrent()->Second();
        return true;
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
        serializer.SetClassDescription( "Triggers if the result changes");
        return serializer;
      }

    }; // template class CriterionResultChangedt< t_ArgumentType, t_ResultType>

    // instantiate s_Instance
    template< typename t_ArgumentType, typename t_ResultType>
    const util::SiPtr< const util::ObjectInterface> CriterionResultChanged< t_ArgumentType, t_ResultType>::s_Instance
    (
      util::Enumerated< CriterionInterface< t_ArgumentType, t_ResultType> >::AddInstance( new CriterionResultChanged< t_ArgumentType, t_ResultType>())
    );

  } // namespace opti
} // namespace bcl

#endif // BCL_OPTI_CRITERION_RESULT_CHANGED_H_
