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

#ifndef BCL_FIND_LOCATOR_CRITERIA_WRAPPER_H_
#define BCL_FIND_LOCATOR_CRITERIA_WRAPPER_H_

// include the namespace header
#include "bcl_find.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_find_locator_criteria_interface.h"
#include "bcl_find_locator_interface.h"
#include "io/bcl_io_serialization.h"
#include "util/bcl_util_enumerated.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace find
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class LocatorCriteriaWrapper
    //! @brief TODO: document
    //! @details TODO: document
    //!
    //! @see @link example_find_locator_criteria_wrapper.cpp @endlink
    //! @author alexanns
    //! @date 01/16/2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_ReturnType, typename t_ArgumentType, typename t_CriteriaType>
    class LocatorCriteriaWrapper :
      public LocatorCriteriaInterface< t_ReturnType, t_ArgumentType, t_CriteriaType>
    {

    private:

    //////////
    // data //
    //////////

      //! Locator without the criteria used in this wrapper
      util::Implementation< LocatorInterface< t_ReturnType, t_ArgumentType> > m_Locator;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      LocatorCriteriaWrapper() :
        m_Locator()
      {
      }

      //! @brief constructor from a locator interface
      LocatorCriteriaWrapper
      (
        const LocatorInterface< t_ReturnType, t_ArgumentType> &LOCATOR
      ) :
        m_Locator( LOCATOR)
      {
      }

      //! @brief constructor from a locator interface
      LocatorCriteriaWrapper
      (
        const util::ShPtr< LocatorInterface< t_ReturnType, t_ArgumentType> > &SP_LOCATOR
      ) :
        m_Locator( *SP_LOCATOR)
      {
      }

      //! @brief virtual copy constructor
      LocatorCriteriaWrapper< t_ReturnType, t_ArgumentType, t_CriteriaType> *Clone() const
      {
        return new LocatorCriteriaWrapper< t_ReturnType, t_ArgumentType, t_CriteriaType>( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief get the name of the object when used in a dynamic context
      //! @return the name of the object when used in a dynamic context
      const std::string &GetAlias() const
      {
        static const std::string s_alias( "LocatorCriteriaWrapper");
        return s_alias;
      }

      //! @brief returns parameters for member data that are set up from the
      //! labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const
      {
        io::Serializer serializer;
        serializer.SetClassDescription( "Locator with criterion.");
        serializer.AddInitializer
        (
          "locator",
          "locator without criterion",
          io::Serialization::GetAgent( &m_Locator)
        );

        return serializer;
      }

      //! locate the t_ReturnType in t_ArgumentType
      //! @param ARGUMENT entity that contains a t_ReturnType
      //! @param CRITERIA
      //! @return returns SiPtr to the locate t_ReturnType
      t_ReturnType Locate( const t_ArgumentType &ARGUMENT, const t_CriteriaType &CRITERIA) const
      {
        // call the locator without the argument
        return m_Locator->Locate( ARGUMENT);
      };

    //////////////////////
    // input and output //
    //////////////////////

    }; // template class LocatorCriteriaWrapper

    // instantiate s_Instance
    template< typename t_ReturnType, typename t_ArgumentType, typename t_CriteriaType>
    const util::SiPtr< const util::ObjectInterface> LocatorCriteriaWrapper< t_ReturnType, t_ArgumentType, t_CriteriaType>::s_Instance
    (
      util::Enumerated< LocatorCriteriaInterface< t_ReturnType, t_ArgumentType, t_CriteriaType> >::AddInstance( new LocatorCriteriaWrapper< t_ReturnType, t_ArgumentType, t_CriteriaType>())
    );

  } // namespace find
} // namespace bcl

#endif //BCL_FIND_LOCATOR_CRITERIA_WRAPPER_H_
