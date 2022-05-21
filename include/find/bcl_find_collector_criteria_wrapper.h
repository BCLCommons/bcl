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

#ifndef BCL_FIND_COLLECTOR_CRITERIA_WRAPPER_H_
#define BCL_FIND_COLLECTOR_CRITERIA_WRAPPER_H_

// include the namespace header
#include "bcl_find.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_find_collector_criteria_interface.h"
#include "bcl_find_collector_interface.h"
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
    //! @class CollectorCriteriaWrapper
    //! @brief Class is used for getting a SiPtrList of substituents (e.x. all CA atoms) from an
    //! argument (e.x. protein model).
    //!
    //! @tparam t_ReturnType is the type of substituents which will be collected
    //! @tparam t_ArgumenType is the type of object from which the substituents will be collected
    //!
    //! @see @link example_find_collector_criteria_wrapper.cpp @endlink
    //! @author alexanns
    //! @date 03/276/2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_ReturnType, typename t_ArgumentType, typename t_CriteriaType>
    class CollectorCriteriaWrapper :
      public CollectorCriteriaInterface< t_ReturnType, t_ArgumentType, t_CriteriaType>
    {

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

      //!< Collector without any criteria to be used in this wrapper
      util::Implementation< CollectorInterface< t_ReturnType, t_ArgumentType> > m_Collector;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! default constructor
      CollectorCriteriaWrapper() :
        m_Collector()
      {
      }

      //! @brief constructor from a collector interface
      CollectorCriteriaWrapper
      (
        const CollectorInterface< t_ReturnType, t_ArgumentType> &COLLECTOR
      ) :
        m_Collector( COLLECTOR)
      {
      }

      //! @brief constructor from a collector interface
      CollectorCriteriaWrapper
      (
        const util::ShPtr< CollectorInterface< t_ReturnType, t_ArgumentType> > &SP_COLLECTOR
      ) :
        m_Collector( *SP_COLLECTOR)
      {
      }

      //! virtual copy constructor
      CollectorCriteriaWrapper< t_ReturnType, t_ArgumentType, t_CriteriaType> *Clone() const
      {
        return new CollectorCriteriaWrapper< t_ReturnType, t_ArgumentType, t_CriteriaType>( *this);
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
        static const std::string s_alias( "CollectorCriteriaWrapper");
        return s_alias;
      }

      //! @brief returns parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const
      {
        io::Serializer serializer;
        serializer.SetClassDescription( "Optimization implementation for Monte Carlo Metropolis algorithms.");
        serializer.AddInitializer
        (
          "collector",
          "collector for substituents",
          io::Serialization::GetAgent( &m_Collector)
        );

        return serializer;
      }

    ////////////////
    // operations //
    ////////////////

      //! Collect the t_ReturnType objects in t_ArgumentType
      //! @param ARGUMENT entity that contains a t_ReturnType
      //! @param CRITERIA
      //! @return returns Group of the collected t_ReturnType objects
      t_ReturnType Collect( const t_ArgumentType &ARGUMENT, const t_CriteriaType &CRITERIA) const
      {
        // call the collector without the argument
        return m_Collector->Collect( ARGUMENT);
      }

    //////////////////////
    // input and output //
    //////////////////////

    }; // template class CollectorCriteriaWrapper

    // instantiate s_Instance
    template< typename t_ReturnType, typename t_ArgumentType, typename t_CriteriaType>
    const util::SiPtr< const util::ObjectInterface> CollectorCriteriaWrapper< t_ReturnType, t_ArgumentType, t_CriteriaType>::s_Instance
    (
      util::Enumerated< CollectorCriteriaInterface< t_ReturnType, t_ArgumentType, t_CriteriaType> >::AddInstance
      (
        new CollectorCriteriaWrapper< t_ReturnType, t_ArgumentType, t_CriteriaType>()
      )
    );

  } // namespace find
} // namespace bcl

#endif //BCL_FIND_COLLECTOR_CRITERIA_WRAPPER_H_
