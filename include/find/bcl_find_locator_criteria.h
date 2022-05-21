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

#ifndef BCL_FIND_LOCATOR_CRITERIA_H_
#define BCL_FIND_LOCATOR_CRITERIA_H_

// include the namespace header
#include "bcl_find.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_find_collector_criteria_interface.h"
#include "bcl_find_locator_criteria_interface.h"
#include "bcl_find_pick_criteria_interface.h"
#include "util/bcl_util_enumerated.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace find
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class LocatorCriteria
    //! @brief template class to combine Collector and Picker templates
    //! @details This generic Locator class combines a Collector
    //!
    //! @see @link example_find_locator_criteria.cpp @endlink
    //! @author linders
    //! @date Aug 23, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_ReturnType, typename t_ArgumentType, typename t_CriteriaType, typename t_IntermediateResultType>
    class LocatorCriteria :
      public LocatorCriteriaInterface< t_ReturnType, t_ArgumentType, t_CriteriaType>
    {

    private:

    //////////
    // data //
    //////////

      //! method for collecting arguments from t_ArgumentType
      util::Implementation< CollectorCriteriaInterface< t_IntermediateResultType, t_ArgumentType, t_CriteriaType> > m_Collector;

      //! method for picking argument from the results of m_Collector
      util::Implementation< PickCriteriaInterface< t_ReturnType, t_IntermediateResultType, t_CriteriaType> > m_Picker;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      LocatorCriteria() :
        m_Collector(),
        m_Picker()
      {
      }

      //! @brief constructor from ShPtrs to a collector and a picker
      //! @param SP_COLLECTOR ShPtr to the collector to be used
      //! @param SP_PICKER ShPtr to the piicker to be used
      LocatorCriteria
      (
        const util::ShPtr< CollectorCriteriaInterface< t_IntermediateResultType, t_ArgumentType, t_CriteriaType> > &SP_COLLECTOR,
        const util::ShPtr< PickCriteriaInterface< t_ReturnType, t_IntermediateResultType, t_CriteriaType> > &SP_PICKER
      ) :
        m_Collector( *SP_COLLECTOR),
        m_Picker( *SP_PICKER)
      {
      }

      //! @brief Clone function
      //! @return pointer to new LocatorCriteria
      LocatorCriteria *Clone() const
      {
        return new LocatorCriteria( *this);
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

      //! @brief returns collector
      //! @return the collector
      const util::Implementation
      <
        CollectorCriteriaInterface< t_IntermediateResultType, t_ArgumentType, t_CriteriaType>
      > &GetCollector() const
      {
        return m_Collector;
      }

      //! @brief get the name of the object when used in a dynamic context
      //! @return the name of the object when used in a dynamic context
      const std::string &GetAlias() const
      {
        static const std::string s_alias( "LocatorCriteria");
        return s_alias;
      }

      //! @brief returns parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const
      {
        io::Serializer serializer;
        serializer.SetClassDescription( "Combine pickers and collectors.");
        serializer.AddInitializer
        (
          "picker",
          "picker for selecting elements",
          io::Serialization::GetAgent( &m_Picker)
        );
        serializer.AddInitializer
        (
          "collector",
          "collector for collecting elements",
          io::Serialization::GetAgent( &m_Collector)
        );

        return serializer;
      }

    ////////////////
    // operations //
    ////////////////

      //! locate the t_ReturnType in t_ArgumentType
      //! @param ARGUMENT entity that contains a t_ReturnType
      //! @param CRITERIA type of criteria
      //! @return returns the located t_ReturnType
      virtual t_ReturnType Locate( const t_ArgumentType &ARGUMENT, const t_CriteriaType &CRITERIA) const
      {
        // collect subset of ARGUMENT that fulfills CRITERIA
        const t_IntermediateResultType collected_subset( m_Collector->Collect( ARGUMENT, CRITERIA));

        // pick one argument from the collected_subset
        const t_ReturnType picked_argument( m_Picker->Pick( collected_subset, CRITERIA));

        // return the picked argument
        return picked_argument;
      }

    //////////////////////
    // input and output //
    //////////////////////

    }; // template class LocatorCriteria

    // instantiate s_Instance
    template
    <
      typename t_ReturnType, typename t_ArgumentType, typename t_CriteriaType, typename t_IntermediateResultType
    >
    const util::SiPtr< const util::ObjectInterface>
    LocatorCriteria< t_ReturnType, t_ArgumentType, t_CriteriaType, t_IntermediateResultType>::s_Instance
    (
      util::Enumerated< LocatorCriteriaInterface< t_ReturnType, t_ArgumentType, t_CriteriaType> >::AddInstance
      (
        new LocatorCriteria< t_ReturnType, t_ArgumentType, t_CriteriaType, t_IntermediateResultType>()
      )
    );

  } // namespace find
} // namespace bcl

#endif // BCL_FIND_LOCATOR_CRITERIA_H_ 
