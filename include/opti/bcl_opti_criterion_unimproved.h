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

#ifndef BCL_OPTI_CRITERION_UNIMPROVED_H_
#define BCL_OPTI_CRITERION_UNIMPROVED_H_

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
    //! @class CriterionUnimproved
    //! @brief CriterionInterface implementation for termination after a certain number of unimproved steps in a row.
    //!
    //! @tparam t_ArgumentType is the type of the approximation argument.
    //! @tparam t_ResultType is the type of the approximation result. A difference has to be supported and a '<' less
    //! than operation
    //!
    //! @see @link example_opti_criterion_unimproved.cpp @endlink
    //! @author fischea
    //! @date Jan 4, 2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_ArgumentType, typename t_ResultType>
    class CriterionUnimproved :
      public CriterionInterface< t_ArgumentType, t_ResultType>
    {

    //////////
    // data //
    //////////

    private:

      //! maximum number of repeats
      size_t m_MaxUnimprovedSteps;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      CriterionUnimproved() :
        m_MaxUnimprovedSteps( 0)
      {
      }

      //! @brief construct from maximum number of repeats and tolerance
      //! @param MAX_NUMBER_UNIMPROVED_STEPS maximum number of unimproved steps
      CriterionUnimproved( const size_t &MAX_NUMBER_UNIMPROVED_STEPS) :
        m_MaxUnimprovedSteps( MAX_NUMBER_UNIMPROVED_STEPS)
      {
      }

      //! @brief Clone function
      //! @return pointer to a new CriterionUnimproved< t_ArgumentType, t_ResultType>
      CriterionUnimproved *Clone() const
      {
        return new CriterionUnimproved< t_ArgumentType, t_ResultType>( *this);
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
        static const std::string s_alias( "ConsecutiveUnimprovedSteps");
        return s_alias;
      }

      //! @brief returns the maximum number of unimproved steps
      //! @return maximum number of unimproved steps
      size_t GetMaxNumberUnimprovedSteps() const
      {
        return m_MaxUnimprovedSteps;
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
        return TRACKER.GetIterationsSinceLastImprovement() > m_MaxUnimprovedSteps;
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
        serializer.SetClassDescription( "Triggers if the result is not improved after the given # of iterations");
        serializer.AddInitializer
        (
          "",
          "",
          io::Serialization::GetAgent( &m_MaxUnimprovedSteps)
        );
        return serializer;
      }

    }; // template class CriterionUnimprovedt< t_ArgumentType, t_ResultType>

    // instantiate s_Instance
    template< typename t_ArgumentType, typename t_ResultType>
    const util::SiPtr< const util::ObjectInterface> CriterionUnimproved< t_ArgumentType, t_ResultType>::s_Instance
    (
      util::Enumerated< CriterionInterface< t_ArgumentType, t_ResultType> >::AddInstance( new CriterionUnimproved< t_ArgumentType, t_ResultType>())
    );

  } // namespace opti
} // namespace bcl

#endif // BCL_OPTI_CRITERION_UNIMPROVED_H_
