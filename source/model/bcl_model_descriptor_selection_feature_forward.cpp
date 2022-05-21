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

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "model/bcl_model_descriptor_selection_feature_forward.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_message.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {
    // add the interface to the set of known implementations
    const util::SiPtr< const util::ObjectInterface> DescriptorSelectionFeatureForward::s_Instance
    (
      util::Enumerated< DescriptorSelectionInterface>::AddInstance
      (
        new DescriptorSelectionFeatureForward()
      )
    );

    //! @brief Clone function
    //! @return pointer to new DescriptorSelectionFeatureForward
    DescriptorSelectionFeatureForward *DescriptorSelectionFeatureForward::Clone() const
    {
      return new DescriptorSelectionFeatureForward();
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &DescriptorSelectionFeatureForward::GetAlias() const
    {
      static const std::string s_Name( "FeatureForwardSelection");
      return s_Name;
    }

    //! @brief assemble all descriptor combinations based on an initial and a total descriptor set as an object label
    //! @param INITIAL initial descriptor set as an object label
    //! @param TOTAL all available descriptor groups in descriptor selection process
    //! @return container with all possible descriptor combinations based on initial descriptor set
    const storage::Vector< util::ObjectDataLabel> DescriptorSelectionFeatureForward::operator()
    (
      const util::ObjectDataLabel &INITIAL,
      const util::ObjectDataLabel &TOTAL
    ) const
    {
      // final vector of descriptor combination pairs
      storage::Vector< util::ObjectDataLabel> qsar_object_list;
      // set for extracting descriptor group already chosen
      storage::Set< util::ObjectDataLabel> initial_set;
      // descriptor groups available for descriptor selection based on initial set
      storage::Vector< util::ObjectDataLabel> entire_minus_initial;

      // iterate over all possible properties
      for
      (
        util::ObjectDataLabel::const_iterator
          itr_prop( INITIAL.Begin()),
          itr_prop_end( INITIAL.End());
        itr_prop != itr_prop_end;
        ++itr_prop
      )
      {
        initial_set.Insert( *itr_prop);
      }

      // iterate over all possible properties
      for
      (
        util::ObjectDataLabel::const_iterator
          itr_prop( TOTAL.Begin()),
          itr_prop_end( TOTAL.End());
        itr_prop != itr_prop_end;
        ++itr_prop
      )
      {
        // only if property is not in initial set it can stay in vector
        if( !initial_set.Contains( *itr_prop))
        {
          entire_minus_initial.PushBack( *itr_prop);
        }
      }

      BCL_MessageDbg
      (
        "entire_minus_initial list of QSAR code: \n" + util::Format()( entire_minus_initial)
      );

      // get the initial arguments
      storage::Vector< util::ObjectDataLabel> initial_labels( INITIAL.GetArguments());

      // iterate over all possible properties
      for
      (
        util::ObjectDataLabel::const_iterator
          itr_prop( entire_minus_initial.Begin()),
          itr_prop_end( entire_minus_initial.End());
        itr_prop != itr_prop_end;
        ++itr_prop
      )
      {
        // add the current label to the initial labels
        initial_labels.PushBack( *itr_prop);

        // make new label containing all the original labels, along with the new one
        util::ObjectDataLabel new_label( INITIAL.GetName(), INITIAL.GetValue(), initial_labels);

        // remove the most recently added label from the initial labels
        initial_labels.PopBack();

        // insert current combination to final list
        qsar_object_list.PushBack( new_label);
      }

      // return list of descriptor combinations
      return qsar_object_list;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer DescriptorSelectionFeatureForward::GetSerializer() const
    {
      return
        io::Serializer().SetCommandLineIdentifier
        (
          "Sequential feature forward selection - starts with single descriptor groups and adds remaining groups to"
          "an initial successor group"
        );
    }

  } // namespace model
} // namespace bcl
