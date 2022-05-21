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
#include "model/bcl_model_descriptor_selection_backward_elimination.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {

    // add the interface to the set of known implementations
    const util::SiPtr< const util::ObjectInterface> DescriptorSelectionBackwardElimination::s_Instance
    (
      util::Enumerated< DescriptorSelectionInterface>::AddInstance
      (
        new DescriptorSelectionBackwardElimination()
      )
    );

    //! @brief Clone function
    //! @return pointer to new DescriptorSelectionBackwardElimination
    DescriptorSelectionBackwardElimination *DescriptorSelectionBackwardElimination::Clone() const
    {
      return new DescriptorSelectionBackwardElimination();
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &DescriptorSelectionBackwardElimination::GetAlias() const
    {
      static const std::string s_Name( "FeatureBackwardElimination");
      return s_Name;
    }

    //! @brief assemble all descriptor combinations based on an initial and a total descriptor set as an object label
    //! @param INITIAL initial descriptor set as an object label
    //! @param TOTAL all available descriptor groups in descriptor selection process
    //! @return container with all possible descriptor combinations based on initial descriptor set
    const storage::Vector< util::ObjectDataLabel> DescriptorSelectionBackwardElimination::operator()
    (
      const util::ObjectDataLabel &INITIAL,
      const util::ObjectDataLabel &TOTAL
    ) const
    {
      // final vector of descriptor combination pairs
      storage::Vector< util::ObjectDataLabel> qsar_object_list;

      // position in qsar_object_list
      size_t position( 0);

      // iterate over all possible properties
      for
      (
        util::ObjectDataLabel::const_iterator itr_prop( INITIAL.Begin()), itr_prop_end( INITIAL.End());
        itr_prop != itr_prop_end;
        ++itr_prop, ++position
      )
      {
        // get the initial arguments
        storage::Vector< util::ObjectDataLabel> initial_labels( INITIAL.Begin(), INITIAL.End());
        initial_labels.RemoveElements( position, 1);

        // make new label containing all the original labels, along with the new one
        util::ObjectDataLabel new_label( INITIAL.GetName(), INITIAL.GetValue(), initial_labels);

        // insert current combination to final list
        qsar_object_list.PushBack( new_label);
      }

      // return list of descriptor combinations
      return qsar_object_list;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer DescriptorSelectionBackwardElimination::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Sequential feature backward elimination - starts with all descriptor groups and removes one group from"
        "an initial successor group"
      );
      return parameters;
    }

  } // namespace model
} // namespace bcl
