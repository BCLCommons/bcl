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
#include "model/bcl_model_retrieve_data_set_encoded_by_model.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"
#include "linal/bcl_linal_vector_reference.h"
#include "model/bcl_model_interface.h"
#include "util/bcl_util_sh_ptr_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {

  //////////
  // data //
  //////////

    // add the interface to the set of known implementations
    const util::SiPtr< const util::ObjectInterface> RetrieveDataSetEncodedByModel::s_Instance
    (
      util::Enumerated< RetrieveDataSetBase>::AddInstance( new RetrieveDataSetEncodedByModel())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief Clone is the virtual Clone constructor
    //! @return a pointer to new LocatorAtomFurthest which is a copy of this
    RetrieveDataSetEncodedByModel *RetrieveDataSetEncodedByModel::Clone() const
    {
      return new RetrieveDataSetEncodedByModel( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &RetrieveDataSetEncodedByModel::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &RetrieveDataSetEncodedByModel::GetAlias() const
    {
      static const std::string s_Name( "EncodeByModel");
      return s_Name;
    }

    //! @brief Set the code / label for the feature (1st part) of the data set
    //! @param CODE the new code
    //! @return the code / label for the feature (1st part) of the data set
    void RetrieveDataSetEncodedByModel::SelectFeatures( const util::ObjectDataLabel &CODE)
    {
      RetrieveDataSetBase::SelectFeatures( CODE);
      m_Features.TryRead( CODE, util::GetLogger());
    }

    //! @brief Set the code / label for the result (2nd part) of the data set
    //! @return the code / label for the result (2nd part) of the data set
    void RetrieveDataSetEncodedByModel::SelectResults( const util::ObjectDataLabel &CODE)
    {
      RetrieveDataSetBase::SelectResults( CODE);
      m_Implementation->SelectResults( CODE);
    }

    //! @brief Set which id columns to retrieve
    //! @param CODE the id column names to retrieve
    void RetrieveDataSetEncodedByModel::SelectIds( const util::ObjectDataLabel &CODE)
    {
      RetrieveDataSetBase::SelectIds( CODE);
      m_Implementation->SelectIds( CODE);
    }

    //! @brief Get the code / label set for the feature (1st part) of the data set with sizes of each property
    //! @return the code / label for the feature (1st part) of the data set with sizes of each property
    //! the feature code set
    FeatureLabelSet RetrieveDataSetEncodedByModel::GetFeatureLabelsWithSizes() const
    {
      storage::Vector< util::ObjectDataLabel> labels;

      const size_t feature_size( m_Model->GetNumberOutputs());

      // create 1 feature label with size 1 for each feature
      storage::Vector< size_t> sizes( feature_size, 1);

      if( GetFeatureCode().GetNumberArguments() != feature_size)
      {
        // the labels will just be the index of each feature
        labels.AllocateMemory( feature_size);
        for( size_t feature_number( 0); feature_number < feature_size; ++feature_number)
        {
          labels.PushBack( util::ObjectDataLabel( util::Format()( feature_number)));
        }
      }
      else
      {
        // use the given labels
        labels = GetFeatureCode().GetArguments();
      }

      return FeatureLabelSet( "Combine", labels, sizes);
    }

    //! @brief Get the code / label for the result (2nd part) of the data set with sizes of each property
    //! @return the code / label for the result (2nd part) of the data set with sizes of each property
    //! the feature code set
    FeatureLabelSet RetrieveDataSetEncodedByModel::GetResultCodeWithSizes() const
    {
      return m_Implementation->GetResultCodeWithSizes();
    }

    //! @brief Get the code / label for the ids of the data set with sizes of each property
    //! @return the code / label for the ids of the data set with sizes of each property
    //! the feature code set
    FeatureLabelSet RetrieveDataSetEncodedByModel::GetIdCodeWithSizes() const
    {
      return m_Implementation->GetIdCodeWithSizes();
    }

    //! @brief get the number of partitions requested by the user, along with the partition ids
    //! @return the number of partitions requested by the user, along with the partition ids
    storage::Pair< size_t, math::RangeSet< size_t> > RetrieveDataSetEncodedByModel::GetNumberPartitionsAndIds() const
    {
      return m_Implementation->GetNumberPartitionsAndIds();
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief generate dataset and apply model to encode given features
    //! @return generated dataset
    util::ShPtr< descriptor::Dataset>
    RetrieveDataSetEncodedByModel::GenerateDataSet()
    {
      util::ShPtr< descriptor::Dataset> dataset_sp( m_Implementation->GenerateDataSet());

      util::ShPtr< FeatureDataSet< float> > features( dataset_sp->GetFeaturesPtr());

      if( !m_Features.GetLabel().IsEmpty())
      {
        linal::Matrix< float> encoded_features( m_Model->operator ()( *features).GetMatrix());
        util::ShPtr< FeatureDataSet< float> > returned_features
        (
          new FeatureDataSet< float>
          (
            encoded_features.GetNumberRows(),
            m_Features.GetOutputFeatureSize(),
            float( 0)
          )
        );
        returned_features->SetFeatureLabelSet( GetFeatureLabelsWithSizes());

        linal::Matrix< float> &returned_features_ref( returned_features->GetRawMatrix());

        for( size_t row_id( 0); row_id < encoded_features.GetNumberRows(); ++row_id)
        {
          m_Features( encoded_features.GetRow( row_id), returned_features_ref.GetRow( row_id).Begin());
        }

        dataset_sp->SetFeatures( returned_features);
      }
      else
      {
        util::ShPtr< FeatureDataSet< float> > encoded_features( m_Model->operator ()( *features).Clone());
        encoded_features->SetFeatureLabelSet( GetFeatureLabelsWithSizes());
        dataset_sp->SetFeatures( encoded_features);
      }

      return dataset_sp;
    } // GenerateDataSet

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer RetrieveDataSetEncodedByModel::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription( "Encode feature of dataset by given model");
      parameters.AddInitializer
      (
        "retriever",
        "dataset retriever to call to get the entire data set",
        io::Serialization::GetAgent( &m_Implementation)
      );

      parameters.AddInitializer
      (
        "storage",
        "encoder storage directory name",
        io::Serialization::GetAgent( &m_ModelRetriever)
      );

      return parameters;
    } // GetParameters

    //! @brief get the nominal (e.g. best estimate without generating the entire dataset) size of the dataset
    //! @return the nominal (e.g. best estimate without generating the entire dataset) size of the dataset
    size_t RetrieveDataSetEncodedByModel::GetNominalSize() const
    {
      return m_Implementation->GetNominalSize();
    }

    //! @brief Set the members of this object from the given LABEL
    //! @param LABEL the label containing members that should be read of this class
    //! @param ERROR_STREAM stream with which to write errors
    bool RetrieveDataSetEncodedByModel::ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM)
    {
      m_Implementation->SelectFeatures( m_ModelRetriever->RetrieveUniqueDescriptor());
      m_Model = m_ModelRetriever->RetrieveEnsemble().FirstElement();

      return true;
    }

  } // namespace model
} // namespace bcl

