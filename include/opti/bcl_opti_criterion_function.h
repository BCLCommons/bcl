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

#ifndef BCL_OPTI_CRITERION_FUNCTION_H_
#define BCL_OPTI_CRITERION_FUNCTION_H_

// include the namespace header
#include "bcl_opti.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_opti_criterion_interface.h"
#include "io/bcl_io_serialization.h"
#include "util/bcl_util_function_interface_serializable.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opti
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class CriterionFunction
    //! @brief CriterionInterface implementation that uses a decider function to evaluate whether the criterion is met
    //!
    //! @tparam t_ArgumentType is the type of the approximation argument
    //! @tparam t_ResultType is the type of the approximation result
    //!
    //! @see @link example_opti_criterion_function.cpp @endlink
    //! @author fischea, mendenjl
    //! @date Sep 03, 2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_ArgumentType, typename t_ResultType>
    class CriterionFunction :
      public CriterionInterface< t_ArgumentType, t_ResultType>
    {

    //////////
    // data //
    //////////

      //! Function to be used to evaluate whether criteria was met
      util::Implementation< util::FunctionInterfaceSerializable< t_ArgumentType, bool> > m_DeciderFunction;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // @brief default constructor
      CriterionFunction() :
        m_DeciderFunction()
      {
      }

      //! @brief construct from maximum number of allowed iterations
      CriterionFunction( const util::FunctionInterfaceSerializable< t_ArgumentType, bool> &DECIIDER_FUNCTION) :
        m_DeciderFunction( DECIIDER_FUNCTION)
      {
      }

      //! @brief Clone function
      //! @return pointer to a new CriterionFunction< t_ArgumentType, t_ResultType>
      CriterionFunction *Clone() const
      {
        return new CriterionFunction< t_ArgumentType, t_ResultType>( *this);
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
        static const std::string s_alias( "Function");
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
        return m_DeciderFunction->operator()( TRACKER.GetCurrent()->First());
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
          "Uses an argument ( " + GetStaticClassName< t_ArgumentType>() + ")-based criterion to trigger"
        );
        serializer.AddInitializer
        (
          "",
          "",
          io::Serialization::GetAgent( &m_DeciderFunction)
        );
        return serializer;
      }

    }; // template class CriterionFunction< t_ArgumentType, t_ResultType>

    // instantiate s_Instance
    template< typename t_ArgumentType, typename t_ResultType>
    const util::SiPtr< const util::ObjectInterface> CriterionFunction< t_ArgumentType, t_ResultType>::s_Instance
    (
      util::Enumerated< CriterionInterface< t_ArgumentType, t_ResultType> >::AddInstance( new CriterionFunction< t_ArgumentType, t_ResultType>())
    );

  } // namespace opti
} // namespace bcl

#endif // BCL_OPTI_CRITERION_FUNCTION_H_
