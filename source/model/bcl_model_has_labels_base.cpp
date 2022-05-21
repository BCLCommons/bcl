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
#include "model/bcl_model_feature_label_set.h"
#include "model/bcl_model_has_labels_base.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {

    //! @brief Set the code / label for the feature (1st part) of the data set
    //! @param CODE the new code
    //! @return the code / label for the feature (1st part) of the data set
    void HasLabelsBase::SetFeatures( const FeatureLabelSet &CODE)
    {
      m_FeatureCodeObject = CODE;
    }

    //! @brief Set the code / label for the result (2nd part) of the data set
    //! @return the code / label for the result (2nd part) of the data set
    void HasLabelsBase::SetResults( const FeatureLabelSet &CODE)
    {
      m_ResultCodeObject = CODE;
    }

    //! @brief Set the code / label for the ids (3rd part) of the data set
    //! @param CODE the new code
    void HasLabelsBase::SetIds( const FeatureLabelSet &CODE)
    {
      m_IdCodeObject = CODE;
    }

    //! @brief Set the code / label for the feature (1st part) of the data set
    //! @return the code / label for the feature (1st part) of the data set
    const FeatureLabelSet &HasLabelsBase::GetFeatureCode() const
    {
      return m_FeatureCodeObject;
    }

    //! @brief Set the code / label for the result (2nd part) of the data set
    //! @return the code / label for the result (2nd part) of the data set
    const FeatureLabelSet &HasLabelsBase::GetResultCode() const
    {
      return m_ResultCodeObject;
    }

    //! @brief Get the code / label for the ids (3rd part) of the data set
    //! @return the code / label for the ids (3rd part) of the data set
    const FeatureLabelSet &HasLabelsBase::GetIdCode() const
    {
      return m_IdCodeObject;
    }

  } // namespace model
} // namespace bcl
