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

#ifndef BCL_MODEL_FEATURE_LABEL_SET_H_
#define BCL_MODEL_FEATURE_LABEL_SET_H_

// include the namespace header
#include "bcl_model.h"

// other forward includes

// includes from bcl - sorted alphabetically
#include "math/bcl_math_range.h"
#include "storage/bcl_storage_map.h"
#include "storage/bcl_storage_pair.h"
#include "storage/bcl_storage_vector.h"
#include "util/bcl_util_implementation.h"
#include "util/bcl_util_implementation_interface.h"
#include "util/bcl_util_object_data_label.h"
#include "util/bcl_util_own_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class FeatureLabelSet
    //! @brief combines the output of SmallMoleculePropertiesInterfaces into one vector
    //!
    //! @see @link example_model_feature_label_set.cpp @endlink
    //! @author mendenjl
    //! @date Apr 13, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API FeatureLabelSet :
      public util::ObjectInterface
    {
    private:

    //////////
    // data //
    //////////

      //! map from object data label to size
      storage::Map< util::ObjectDataLabel, math::Range< size_t> > m_PropertiesToRanges;

      //! All sub-labels in the order in which they were defined
      storage::Vector< util::ObjectDataLabel> m_OrderedProperties;

      //! Sizes of each property
      storage::Vector< size_t> m_PropertiesSizes;

      //! the size of the combined properties
      size_t m_Size;

      //! Primary name for the label set; usually is a property that combines multiple descriptors, like Combine
      std::string m_OuterName;

      //! Namespace of all the given labels
      util::Implementation< util::ImplementationInterface> m_ImplementationInterface;

      //! shared pointer to the fully split feature label set; necessary if the original feature label set had
      //! partial descriptors
      mutable util::OwnPtr< FeatureLabelSet> m_SplitFeatures;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      FeatureLabelSet();

      //! @brief constructor from outer name
      //! @param NAME name for the set, usually a descriptor that takes a list of other descriptors
      //! @param IMPL implementation interface implementation
      FeatureLabelSet
      (
        const std::string &NAME,
        const util::Implementation< util::ImplementationInterface> &IMPL
      );

      //! @brief constructor from members
      //! @param NAME name for the set, usually a descriptor that takes a list of other descriptors
      //! @param PROPERTIES properties of the feature
      //! @param PROPERTY_SIZES sizes of those properties
      //! @param IMPL implementation interface implementation
      FeatureLabelSet
      (
        const std::string &NAME,
        const storage::Vector< util::ObjectDataLabel> &PROPERTIES,
        const storage::Vector< size_t> &PROPERTY_SIZES,
        const util::Implementation< util::ImplementationInterface> &IMPL =
          util::Implementation< util::ImplementationInterface>()
      );

      //! @brief virtual copy constructor
      FeatureLabelSet *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief return the alias / name used when constructing this object
      //! @returnthe alias / name used when constructing this object
      const std::string &GetAlias() const;

      //! @brief return a label with the descriptors in this
      //! @return a label with the descriptors in this
      util::ObjectDataLabel GetLabel() const;

      //! @brief return the data label
      //! @return data label as string
      std::string GetString() const;

      //! @brief push back a property and its expected size
      //! @param PROPERTY_LABEL label describing the property
      //! @param PROPERTY_SIZE size of the label
      void PushBack
      (
        const util::ObjectDataLabel &PROPERTY_LABEL,
        const size_t &PROPERTY_SIZE
      );

      //! @param LABEL label of the desired property
      //! @return segment(s) of feature vector used by the given property
      storage::Vector< size_t> GetPropertyIndices( const util::ObjectDataLabel &LABEL) const;

      //! @brief Get a vector containing all member data from this property as separate data labels
      const storage::Vector< util::ObjectDataLabel> &GetMemberLabels() const;

      //! @brief get the sizes of each property
      //! @return sizes of each property, in the same order as GetMemberLabels
      const storage::Vector< size_t> &GetPropertySizes() const;

      //! @return number of values represented by this feature label
      size_t GetSize() const
      {
        return m_Size;
      }

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief create a feature label subset
      //! @param SUB_FEATURES indices of the complete feature to keep
      //! @return feature label set with the given sub features
      //! @note: if a descriptor has multiple return values and must be split, it is assumed that there is a
      //! meta-descriptor named Partial, that takes another descriptor and a parameter called indices, which contains
      //! the indices of the feature labels to keep.
      FeatureLabelSet CreateSubFeatureLabelSet( const storage::Vector< size_t> &SUB_FEATURES) const;

      //! @brief create a feature label subset
      //! @param SKIP_ZERO_LENGTH_DESCRIPTORS true to ignore zero-length descriptors during the split. In this case, all
      //!        descriptors in the returned labels will have size == 1
      //! @return feature label set with each multi-column descriptor split using the Partial descriptor (vector index)
      //! @note: if a descriptor has multiple return values and must be split, it is assumed that there is a
      //! meta-descriptor named Partial, that takes another descriptor and a parameter called indices, which contains
      //! the indices of the feature labels to keep.
      FeatureLabelSet SplitFeatureLabelSet( const bool &SKIP_ZERO_LENGTH_DESCRIPTORS = false) const;

      //! @brief get common overlap of features in this and given FeatureLabelSets
      FeatureLabelSet GetCommonFeatures( const FeatureLabelSet &COMPARE) const;

      //! @brief merge two object data labels, respecting split descriptors
      //! @param LABEL_A, LABEL_B the object data labels of interest
      //! @return the merged labels
      static util::ObjectDataLabel MergeConsideringPartials
      (
        const util::ObjectDataLabel &LABEL_A,
        const util::ObjectDataLabel &LABEL_B
      );

      //! @brief helper function to decompose a partial into the innermost object data label and the vector of indices
      //! @param LABEL the partial label of interest
      //! @return a pair containing the main descriptors object data label and the indices for the partial
      static storage::Pair< util::ObjectDataLabel, storage::Vector< size_t> >
        DecomposePartial( const util::ObjectDataLabel &LABEL);

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

    protected:

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief get a vector with all known synonyms for a given label, created from trying to update the data label
      //!        with all template instances of descriptor::Base
      //! @param LABEL the label to retrieve all synonyms for. These are primarily necessary when
      static storage::Set< util::ObjectDataLabel> GetSynonyms( const util::ObjectDataLabel &LABEL);

      //! @brief clean a label that might have multiple levels of partials
      //! @param LABEL the old label to consider
      static util::ObjectDataLabel CollapsePartials( const util::ObjectDataLabel &LABEL);

      //! @brief retrieve the indices of a sub-property
      //! @param LABEL property of the form Partial(x,indices(y1,y2,...,yn))
      //! @return y1,y2,...yn of the property x, if it could be located
      storage::Vector< size_t> GetSubPropertyIndices( const util::ObjectDataLabel &LABEL) const;

    }; // class FeatureLabelSet

  } // namespace model
} // namespace bcl

#endif // BCL_MODEL_FEATURE_LABEL_SET_H_
