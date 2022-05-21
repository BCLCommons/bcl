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

#ifndef BCL_CHEMISTRY_SMALL_MOLECULE_QSAR_STORAGE_FILE_H_
#define BCL_CHEMISTRY_SMALL_MOLECULE_QSAR_STORAGE_FILE_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically
#include "descriptor/bcl_descriptor.fwd.hh"
#include "linal/bcl_linal.fwd.hh"
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_conformation_interface.h"
#include "bcl_chemistry_descriptor_dimension.h"
#include "bcl_chemistry_molecule_storage_file.h"
#include "descriptor/bcl_descriptor_dataset_builder.h"
#include "math/bcl_math_range_set.h"
#include "model/bcl_model_retrieve_data_set_base.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SmallMoleculeQsarStorageFile
    //! @brief is a handler for storing SmallMolecule properties used as molecular descriptors and
    //!        biological activity data
    //!
    //! In an file the key corresponds to the index of molecules which is auto incremented and is only numeric
    //! when initialized as attached, the largest key is located
    //! when initialized as created or overwrite, the key starts with 0
    //!
    //! @see @link example_chemistry_small_molecule_qsar_storage_file.cpp @endlink
    //! @author butkiem1
    //! @date Aug 19, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SmallMoleculeQsarStorageFile :
      public MoleculeStorageFile,
      public model::RetrieveDataSetBase
    {

    private:

    //////////
    // data //
    //////////

      storage::Map< size_t, size_t> m_FeatureStartIdToKey; //!< logical feature start id to key in model storage
      //! Method by which to determine the dimensionality / type of descriptors generated
      DescriptorDimensionEnum m_DescriptorType;
      mutable util::ShPtr< descriptor::DatasetBuilder< AtomConformationalInterface> > m_Builder;

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      SmallMoleculeQsarStorageFile();

      //! @brief copy constructor
      SmallMoleculeQsarStorageFile( const SmallMoleculeQsarStorageFile &PARENT);

      //! @brief Clone function
      //! @return pointer to new SmallMoleculeQsarStorageFile
      SmallMoleculeQsarStorageFile *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief returns the name used for this class in an object data label
      //! @return the name used for this class in an object data label
      const std::string &GetAlias() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief Get the code / label set for the feature (1st part) of the data set with sizes of each property
      //! @return the code / label for the feature (1st part) of the data set with sizes of each property
      //! the feature code set
      model::FeatureLabelSet GetFeatureLabelsWithSizes() const;

      //! @brief Get the code / label for the result (2nd part) of the data set with sizes of each property
      //! @return the code / label for the result (2nd part) of the data set with sizes of each property
      //! the feature code set
      model::FeatureLabelSet GetResultCodeWithSizes() const;

      //! @brief Get the code / label for the ids of the data set with sizes of each property
      //! @return the code / label for the ids of the data set with sizes of each property
      //! the feature code set
      model::FeatureLabelSet GetIdCodeWithSizes() const;

      //! @brief Set the code / label for the feature (1st part) of the data set
      //! @param CODE the new code
      void SelectFeatures( const util::ObjectDataLabel &CODE);

      //! @brief Set the code / label for the result (2nd part) of the data set
      //! @param CODE the new code
      void SelectResults( const util::ObjectDataLabel &CODE);

      //! @brief Set the code / label for the ids (3rd part) of the data set
      //! @param CODE the new code
      void SelectIds( const util::ObjectDataLabel &CODE);

      //! @brief get whether dataset generation requires labels
      //! @return true if dataset generation requires labels
      bool RequiresFeatureLabels() const
      {
        return true;
      }

      //! @brief get whether dataset generation requires result labels
      //! @return true if dataset generation requires result labels
      bool RequiresResultLabels() const
      {
        return true;
      }

      //! @brief generate dataset from a set of ranges
      //! @return generated dataset
      util::ShPtr< descriptor::Dataset> GenerateDataSet();

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

      //! @brief get the number of partitions requested by the user, along with the partition ids
      //! @return the number of partitions requested by the user, along with the partition ids
      storage::Pair< size_t, math::RangeSet< size_t> > GetNumberPartitionsAndIds() const;

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

      //! @brief get the nominal (e.g. best estimate without generating the entire dataset) size of the dataset
      //! @return the nominal (e.g. best estimate without generating the entire dataset) size of the dataset
      size_t GetNominalSize() const;

      //! @brief test whether this retriever can generate sub-ranges of datasets without loading the entire dataset
      //! @return true if this retriever can generate sub-ranges of datasets without loading the entire dataset
      bool SupportsEfficientSubsetLoading() const
      {
        return true;
      }

      //! @brief load a range of data from the dataset
      //! @param SUBSET the range of data to load
      //! @param FEATURES_STORAGE where to store features that are loaded, must be large enough to hold the subset without resizing
      //! @param RESULTS_STORAGE where to store the corresponding results, must be large enough to hold the subset without resizing
      //! @param START_FEATURE_NUMBER position to store the first feature in FEATURES_STORAGE
      //! @return # of features actually loaded
      //! Note: Implementations should overload this and SupportsEfficientSubsetLoading together
      size_t GenerateDataSubset
      (
        const math::Range< size_t> &SUBSET,
        linal::MatrixInterface< float> &FEATURES_STORAGE,
        linal::MatrixInterface< float> &RESULTS_STORAGE,
        linal::MatrixInterface< char> &IDS_STORAGE,
        const size_t &START_FEATURE_NUMBER
      );

      //! @brief load a range of data from the dataset, given a particular dataset builder
      //! @param SUBSET the range of data to load
      //! @param FEATURES_STORAGE where to store features that are loaded, must be large enough to hold the subset without resizing
      //! @param RESULTS_STORAGE where to store the corresponding results, must be large enough to hold the subset without resizing
      //! @param START_FEATURE_NUMBER position to store the first feature in FEATURES_STORAGE
      //! @param BUILDER the dataset builder to use
      //! @return # of features actually loaded
      size_t GenerateDataSubsetGivenBuilder
      (
        const math::Range< size_t> &SUBSET,
        linal::MatrixInterface< float> &FEATURES_STORAGE,
        linal::MatrixInterface< float> &RESULTS_STORAGE,
        linal::MatrixInterface< char> &IDS_STORAGE,
        const size_t &START_FEATURE_NUMBER,
        descriptor::DatasetBuilder< AtomConformationalInterface> &BUILDER
      );

      //! @brief fix the descriptor type if it was specified over the command line
      void FixDescriptorType();

    }; // class SmallMoleculeQsarStorageFile

  } // namespace chemistry

} // namespace bcl

#endif // BCL_CHEMISTRY_SMALL_MOLECULE_QSAR_STORAGE_FILE_H_
