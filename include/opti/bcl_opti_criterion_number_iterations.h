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

#ifndef BCL_OPTI_CRITERION_NUMBER_ITERATIONS_H_
#define BCL_OPTI_CRITERION_NUMBER_ITERATIONS_H_

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
    //! @class CriterionNumberIterations
    //! @brief CriterionInterface implementation for terminating after a specified number of iterations
    //!
    //! @tparam t_ArgumentType is the type of the approximation argument
    //! @tparam t_ResultType is the type of the approximation result
    //!
    //! @see @link example_opti_criterion_number_iterations.cpp @endlink
    //! @author fischea
    //! @date Dec 12, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_ArgumentType, typename t_ResultType>
    class CriterionNumberIterations :
      public CriterionInterface< t_ArgumentType, t_ResultType>
    {

    //////////
    // data //
    //////////

    private:

      //! maximum number of allowed iterations
      size_t m_MaxNumberIterations;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // @brief default constructor
      CriterionNumberIterations() :
        m_MaxNumberIterations( 0)
      {
      }

      //! @brief construct from maximum number of allowed iterations
      //! @param MAX_NUMBER_ITERATIONS maximum number of allowed iterations
      CriterionNumberIterations( const size_t &MAX_NUMBER_ITERATIONS) :
        m_MaxNumberIterations( MAX_NUMBER_ITERATIONS)
      {
      }

      //! @brief Clone function
      //! @return pointer to a new CriterionNumberIterations< t_ArgumentType, t_ResultType>
      CriterionNumberIterations *Clone() const
      {
        return new CriterionNumberIterations< t_ArgumentType, t_ResultType>( *this);
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
        static const std::string s_alias( "Iterations");
        return s_alias;
      }

      //! @brief returns the maximum number of allowed iterations
      //! @return maximum number of allowed iterations
      size_t GetMaxNumberIterations() const
      {
        return m_MaxNumberIterations;
      }

      //! @brief sets the maximum number of allowed iterations
      //! @param MAX_NUMBER_ITERATIONS maximum number of allowed iterations
      void SetMaxNumberIterations( const size_t &MAX_NUMBER_ITERATIONS)
      {
        m_MaxNumberIterations = MAX_NUMBER_ITERATIONS;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief returns whether the maximum number of iterations has been observed for the given tracker
      //! @param TRACKER tracker to evaluate the criterion for
      //! @return true, if any of the termination criteria are met yet
      bool CriteriaMet( const Tracker< t_ArgumentType, t_ResultType> &TRACKER) const
      {
        return TRACKER.GetIteration() >= m_MaxNumberIterations;
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
        serializer.SetClassDescription( "Triggers after the given number of iterations");
        serializer.AddInitializer
        (
          "",
          "",
          io::Serialization::GetAgent( &m_MaxNumberIterations)
        );
        return serializer;
      }

    }; // template class CriterionNumberIterations< t_ArgumentType, t_ResultType>

    // instantiate s_Instance
    template< typename t_ArgumentType, typename t_ResultType>
    const util::SiPtr< const util::ObjectInterface> CriterionNumberIterations< t_ArgumentType, t_ResultType>::s_Instance
    (
      util::Enumerated< CriterionInterface< t_ArgumentType, t_ResultType> >::AddInstance( new CriterionNumberIterations< t_ArgumentType, t_ResultType>())
    );

  } // namespace opti
} // namespace bcl

#endif // BCL_OPTI_CRITERION_NUMBER_ITERATIONS_H_
