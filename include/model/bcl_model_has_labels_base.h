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

#ifndef BCL_MODEL_HAS_LABELS_BASE_H_
#define BCL_MODEL_HAS_LABELS_BASE_H_

// include the namespace header
#include "bcl_model.h"
#include "bcl_model_feature_label_set.h"
#include "util/bcl_util_object_data_label.h"
#include "util/bcl_util_object_interface.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class HasLabelsBase
    //! @brief This class provides basic functionality for a approximation scheme to allow access to data set
    //!        related object data labels
    //!
    //! @remarks example unnecessary
    //! @author butkiem1, mendenjl
    //! @date Dec 04, 2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API HasLabelsBase :
      public virtual util::ObjectInterface
    {

    //////////
    // data //
    //////////

    protected:

      //! id code object
      FeatureLabelSet m_IdCodeObject;
      //! feature code object
      FeatureLabelSet m_FeatureCodeObject;
      //! result code object
      FeatureLabelSet m_ResultCodeObject;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

    public:

      //! @brief Set the code / label for the feature (1st part) of the data set
      //! @param CODE the new code
      virtual void SetFeatures( const FeatureLabelSet &CODE);

      //! @brief Set the code / label for the result (2nd part) of the data set
      //! @param CODE the new code
      virtual void SetResults( const FeatureLabelSet &CODE);

      //! @brief Set the code / label for the ids (3rd part) of the data set
      //! @param CODE the new code
      virtual void SetIds( const FeatureLabelSet &CODE);

      //! @brief Get the code / label for the feature (1st part) of the data set
      //! @return the code / label for the feature (1st part) of the data set
      virtual const FeatureLabelSet &GetFeatureCode() const;

      //! @brief Get the code / label for the result (2nd part) of the data set
      //! @return the code / label for the result (2nd part) of the data set
      virtual const FeatureLabelSet &GetResultCode() const;

      //! @brief Get the code / label for the ids (3rd part) of the data set
      //! @return the code / label for the ids (3rd part) of the data set
      virtual const FeatureLabelSet &GetIdCode() const;

    }; // class HasLabelsBase

  } // namespace model
} // namespace bcl

#endif // BCL_MODEL_HAS_LABELS_BASE_H_
