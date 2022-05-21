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

#ifndef BCL_OPTI_CRITERION_DIVERGENCE_ARGUMENT_H_
#define BCL_OPTI_CRITERION_DIVERGENCE_ARGUMENT_H_

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
    //! @class CriterionDivergenceArgument
    //! @brief CriterionInterface implementation that terminates after a repetitive difference between previous and
    //! current argument larger that given interval width
    //!
    //! @tparam t_ArgumentType is the type of the approximation argument. a difference has to be supported and a '<'
    //! less than operation
    //! @tparam t_ResultType is the type of the approximation result.
    //!
    //! @see @link example_opti_criterion_divergence_argument.cpp @endlink
    //! @author woetzen, fischea
    //! @date Dec 13, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_ArgumentType, typename t_ResultType>
    class CriterionDivergenceArgument :
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

      //! tolerance to consider arguments as repeats
      t_ArgumentType m_IntervalWidth;

      //! previous argument result pair
      mutable t_ArgumentType m_PreviousArgument;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      CriterionDivergenceArgument() :
        m_NumberRepeats( util::GetUndefined< size_t>()),
        m_MaxNumberRepeats( 0),
        m_IntervalWidth( 0.0)
      {
      }

      //! @brief construct from maximum number of repeats and tolerance
      //! @param MAX_NUMBER_REPEATS maximum number of repeats
      //! @param TOLERANCE tolerance to consider arguments as repeats
      CriterionDivergenceArgument
      (
        const size_t &MAX_NUMBER_REPEATS,
        const t_ArgumentType &INTERVAL_WIDTH
      ) :
        m_NumberRepeats( util::GetUndefined< size_t>()),
        m_MaxNumberRepeats( MAX_NUMBER_REPEATS),
        m_IntervalWidth( INTERVAL_WIDTH)
      {
      }

      //! @brief Clone function
      //! @return pointer to a new CriterionDivergenceArgument< t_ArgumentType, t_ResultType>
      CriterionDivergenceArgument *Clone() const
      {
        return new CriterionDivergenceArgument< t_ArgumentType, t_ResultType>( *this);
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
        static const std::string s_alias( "DivergenceArgument");
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
      t_ArgumentType GetIntervalWidth() const
      {
        return m_IntervalWidth;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief returns true if maximum number of iterations with difference of two successive results larger than the
      //! interval width has been observed
      //! @param TRACKER tracker to evaluate the criterion for
      //! @return true, true if maximum number of iterations with difference of two successive results larger than the
      //! interval width has been observed
      bool CriteriaMet( const Tracker< t_ArgumentType, t_ResultType> &TRACKER) const
      {
        // get current argument result pair from the tracker
        const util::ShPtr< storage::Pair< t_ArgumentType, t_ResultType> > &sp_current( TRACKER.GetCurrent());

        if( !sp_current.IsDefined())
        {
          return false;
        }

        // if no previous argument result pair has been observed, the criterion is not met
        if( !util::IsDefined( m_NumberRepeats))
        {
          m_PreviousArgument = sp_current->First();
          m_NumberRepeats = 0;
          return false;
        }

        // check whether the new argument lies within the tolerance
        if( math::Absolute( sp_current->First() - m_PreviousArgument) > m_IntervalWidth)
        {
          ++m_NumberRepeats;
        }
        else
        {
          m_NumberRepeats = 0;
        }

        // set previous tracker result to current item
        m_PreviousArgument = sp_current->First();

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
          "Criteria that the argument diverges its previous value for a given number of iterations in a row"
        );
        serializer.AddInitializer
        (
          "repeats",
          "Required # of consecutive iterations in which the argument has changed more than the tolerance",
          io::Serialization::GetAgent( &m_MaxNumberRepeats)
        );
        serializer.AddInitializer
        (
          "tolerance",
          "tolerance under which arguments are considered as repeats",
          io::Serialization::GetAgent( &m_IntervalWidth)
        );
        serializer.AddDataMember( "current repetitions", io::Serialization::GetAgent( &m_NumberRepeats));
        serializer.AddDataMember( "last argument", io::Serialization::GetAgent( &m_PreviousArgument));
        return serializer;
      }

    }; // template class CriterionDivergenceArgument< t_ArgumentType, t_ResultType>

    // instantiate s_Instance
    template< typename t_ArgumentType, typename t_ResultType>
    const util::SiPtr< const util::ObjectInterface> CriterionDivergenceArgument< t_ArgumentType, t_ResultType>::s_Instance
    (
      util::Enumerated< CriterionInterface< t_ArgumentType, t_ResultType> >::AddInstance( new CriterionDivergenceArgument< t_ArgumentType, t_ResultType>())
    );

  } // namespace opti
} // namespace bcl

#endif // BCL_OPTI_CRITERION_DIVERGENCE_ARGUMENT_H_
