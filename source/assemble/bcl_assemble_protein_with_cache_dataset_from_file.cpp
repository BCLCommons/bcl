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
#include "assemble/bcl_assemble_protein_with_cache_dataset_from_file.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_classes.h"
#include "descriptor/bcl_descriptor_dataset_builder.h"
#include "io/bcl_io_serialization.h"
#include "pdb/bcl_pdb_factory.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {

  //////////
  // data //
  //////////

    const util::SiPtr< const util::ObjectInterface> ProteinWithCacheDatasetFromFile::s_SequenceInstance
    (
      util::Enumerated< model::RetrieveDataSetBase>::AddInstance( new ProteinWithCacheDatasetFromFile( false))
    );
    const util::SiPtr< const util::ObjectInterface> ProteinWithCacheDatasetFromFile::s_ProteinInstance
    (
      util::Enumerated< model::RetrieveDataSetBase>::AddInstance( new ProteinWithCacheDatasetFromFile( true))
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor, takes the desired AA class and whether or not to require coordinates
    ProteinWithCacheDatasetFromFile::ProteinWithCacheDatasetFromFile( const bool &REQUIRE_COORDINATES) :
      ProteinWithCacheStorageFile( REQUIRE_COORDINATES),
      m_Builder( new descriptor::DatasetBuilder< biol::AABase>())
    {
    }

    //! @brief copy constructor, clones builder
    ProteinWithCacheDatasetFromFile::ProteinWithCacheDatasetFromFile( const ProteinWithCacheDatasetFromFile &PARENT) :
      ProteinWithCacheStorageFile( PARENT),
      model::RetrieveDataSetBase( PARENT),
      m_FeatureStartIdToKey( PARENT.m_FeatureStartIdToKey),
      m_Builder
      (
        new descriptor::DatasetBuilder< biol::AABase>
        (
          PARENT.GetFeatureCode(),
          PARENT.GetResultCode(),
          PARENT.GetIdCode()
        )
      )
    {
    }

    //! @brief Clone function
    //! @return pointer to new ProteinWithCacheDatasetFromFile
    ProteinWithCacheDatasetFromFile *ProteinWithCacheDatasetFromFile::Clone() const
    {
      return new ProteinWithCacheDatasetFromFile( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ProteinWithCacheDatasetFromFile::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &ProteinWithCacheDatasetFromFile::GetAlias() const
    {
      return ProteinWithCacheStorageFile::GetAlias();
    }

    //! @brief Set the code / label for the feature (1st part) of the data set
    //! @param CODE the new code
    void ProteinWithCacheDatasetFromFile::SelectFeatures( const util::ObjectDataLabel &CODE)
    {
      model::RetrieveDataSetBase::SelectFeatures( CODE);
      // this is not strictly necessary, but it allows a user to ask for help over the command line for this retriever
      util::Implementation< descriptor::Base< biol::AABase, float> > property( GetFeatureCode());
      m_Builder->SetFeatureCode( GetFeatureCode());
    }

    //! @brief Set the code / label for the result (2nd part) of the data set
    //! @param CODE the new code
    void ProteinWithCacheDatasetFromFile::SelectResults( const util::ObjectDataLabel &CODE)
    {
      model::RetrieveDataSetBase::SelectResults( CODE);
      m_FeatureStartIdToKey.Reset();
      m_Builder->SetResultsCode( GetResultCode());
      const descriptor::Type type( m_Builder->GetType());
      if( type.GetDimension() == size_t( 0))
      {
        // add up the size of each, convoluted by the type
        size_t key_number( 0);
        for( const size_t n_proteins( ProteinWithCacheStorageFile::GetSize()); key_number < n_proteins; ++key_number)
        {
          m_FeatureStartIdToKey[ key_number] = key_number;
        }
        m_FeatureStartIdToKey[ key_number] = key_number;
      }
      else
      {
        // retrieve the size of each key
        storage::Vector< std::string> all_keys( ProteinWithCacheStorageFile::GetAllKeys());
        // add up the size of each, convoluted by the type
        size_t total_features( 0), key_number( 0);
        util::GetLogger() << "Computing total # of features\n";
        for
        (
          storage::Vector< std::string>::const_iterator itr( all_keys.Begin()), itr_end( all_keys.End());
          itr != itr_end;
          ++itr, ++key_number
        )
        {
          const size_t n_new_features( type.GetNumberFeatures( ProteinWithCacheStorageFile::GetKeySize( *itr)));
          if( n_new_features)
          {
            m_FeatureStartIdToKey[ total_features] = key_number;
            total_features += n_new_features;
            util::GetLogger().LogStatus
            (
              "File #" + util::Format()( key_number) + " / " + util::Format()( all_keys.GetSize())
              + " (" + *itr + "), # features: " + util::Format()( n_new_features)
              + "; total features so far: " + util::Format()( total_features)
            );
          }
        }
        m_FeatureStartIdToKey[ total_features] = key_number;
      }
    }

    //! @brief Set the code / label for the ids (3rd part) of the data set
    //! @param CODE the new code
    void ProteinWithCacheDatasetFromFile::SelectIds( const util::ObjectDataLabel &CODE)
    {
      model::RetrieveDataSetBase::SelectIds( CODE);
      // this is not strictly necessary, but it allows a user to ask for help over the command line for this retriever
      util::Implementation< descriptor::Base< biol::AABase, char> > property( GetIdCode());
      m_Builder->SetIdCode( CODE);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief Get the code / label set for the feature (1st part) of the data set with sizes of each property
    //! @return the code / label for the feature (1st part) of the data set with sizes of each property
    //! the feature code set
    model::FeatureLabelSet ProteinWithCacheDatasetFromFile::GetFeatureLabelsWithSizes() const
    {
      return m_Builder.IsDefined() ? m_Builder->GetFeatureCode().GetLabelsWithSizes() : model::FeatureLabelSet();
    }

    //! @brief Get the code / label for the result (2nd part) of the data set with sizes of each property
    //! @return the code / label for the result (2nd part) of the data set with sizes of each property
    //! the feature code set
    model::FeatureLabelSet ProteinWithCacheDatasetFromFile::GetResultCodeWithSizes() const
    {
      return m_Builder.IsDefined() ? m_Builder->GetResultCode().GetLabelsWithSizes() : model::FeatureLabelSet();
    }

    //! @brief Get the code / label for the ids of the data set with sizes of each property
    //! @return the code / label for the ids of the data set with sizes of each property
    //! the feature code set
    model::FeatureLabelSet ProteinWithCacheDatasetFromFile::GetIdCodeWithSizes() const
    {
      return m_Builder.IsDefined() ? m_Builder->GetIdCode().GetLabelsWithSizes() : model::FeatureLabelSet();
    }

    //! @brief generate dataset from a set of ranges
    //! @return generated dataset
    util::ShPtr< descriptor::Dataset> ProteinWithCacheDatasetFromFile::GenerateDataSet()
    {
      pdb::Factory::GetFlagConvertToNaturalAAType()->SetFlag();

      // determine the total size of the data set
      const size_t total_dataset_size( m_FeatureStartIdToKey.ReverseBegin()->first);

      if( !m_Builder.IsDefined())
      {
        m_Builder =
          util::ShPtr< descriptor::DatasetBuilder< biol::AABase> >
          (
            new descriptor::DatasetBuilder< biol::AABase>( GetFeatureCode(), GetResultCode(), GetIdCode())
          );
      }

      // dataset to hold all results
      const size_t nominal_size( GetNominalSize());
      util::ShPtr< descriptor::Dataset> dataset_sp
      (
        new descriptor::Dataset
        (
          nominal_size,
          GetFeatureLabelsWithSizes(),
          GetResultCodeWithSizes(),
          GetIdCodeWithSizes()
        )
      );
      linal::MatrixReference< float> features( dataset_sp->GetFeaturesReference());
      linal::MatrixReference< float> results( dataset_sp->GetResultsReference());
      linal::MatrixReference< char> ids( dataset_sp->GetIdsReference());

      const size_t nr_features
      (
        GenerateDataSubsetGivenBuilder
        (
          math::Range< size_t>( 0, total_dataset_size),
          features,
          results,
          ids,
          0,
          *m_Builder
        )
      );

      // remove unused rows
      dataset_sp->ShrinkRows( nr_features);

      // return the generated dataset
      return dataset_sp;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ProteinWithCacheDatasetFromFile::Read( std::istream &ISTREAM)
    {
      // read member
      ProteinWithCacheStorageFile::Read( ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &ProteinWithCacheDatasetFromFile::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member
      ProteinWithCacheStorageFile::Write( OSTREAM, INDENT) << '\n';

      // end
      return OSTREAM;
    }

    //! @brief get the number of partitions requested by the user, along with the partition ids
    //! @return the number of partitions requested by the user, along with the partition ids
    storage::Pair< size_t, math::RangeSet< size_t> > ProteinWithCacheDatasetFromFile::GetNumberPartitionsAndIds() const
    {
      return storage::Pair< size_t, math::RangeSet< size_t> >( 1, math::RangeSet< size_t>( math::Range< size_t>( 0, 0)));
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer ProteinWithCacheDatasetFromFile::GetSerializer() const
    {
      io::Serializer member_data( ProteinWithCacheStorageFile::GetSerializer());
      member_data.SetClassDescription
      (
        (
          ProteinWithCacheStorageFile::m_RequireCoordinates
          ? "Calculates descriptors from PDB entries (as sequences of " + ProteinWithCacheStorageFile::m_Class.GetName() + ")"
          : "Retrieves sequences of " + ProteinWithCacheStorageFile::m_Class.GetName() + " from PDBs or FASTAs"
        ) +
        " from a directory or subdirectory"
      );

      return member_data;
    } // GetParameters

    //! @brief get the nominal (e.g. best estimate without generating the entire dataset) size of the dataset
    //! @return the nominal (e.g. best estimate without generating the entire dataset) size of the dataset
    size_t ProteinWithCacheDatasetFromFile::GetNominalSize() const
    {
      // get the last value out of the map, if there are any values
      if( m_FeatureStartIdToKey.IsEmpty())
      {
        BCL_MessageCrt( "Warning, called GetNominalSize before calling SelectResults");
        return 0;
      }

      return m_FeatureStartIdToKey.ReverseBegin()->first;
    }

    //! @brief load a range of data from the dataset
    //! @param SUBSET the range of data to load
    //! @param FEATURES_STORAGE where to store features that are loaded, must be large enough to hold the subset without resizing
    //! @param RESULTS_STORAGE where to store the corresponding results, must be large enough to hold the subset without resizing
    //! @param START_FEATURE_NUMBER position to store the first feature in FEATURES_STORAGE
    //! @return # of features actually loaded
    //! Note: Implementations should overload this and SupportsEfficientSubsetLoading together
    size_t ProteinWithCacheDatasetFromFile::GenerateDataSubset
    (
      const math::Range< size_t> &SUBSET,
      linal::MatrixInterface< float> &FEATURES_STORAGE,
      linal::MatrixInterface< float> &RESULTS_STORAGE,
      linal::MatrixInterface< char> &IDS_STORAGE,
      const size_t &START_FEATURE_NUMBER
    )
    {
      pdb::Factory::GetFlagConvertToNaturalAAType()->SetFlag();

      if( !m_Builder.IsDefined())
      {
        m_Builder =
          util::ShPtr< descriptor::DatasetBuilder< biol::AABase> >
          (
            new descriptor::DatasetBuilder< biol::AABase>( GetFeatureCode(), GetResultCode(), GetIdCode())
          );
      }

      // generate a dataset of protein codes from the ranges given in the range set
      return
        GenerateDataSubsetGivenBuilder
        (
          SUBSET,
          FEATURES_STORAGE,
          RESULTS_STORAGE,
          IDS_STORAGE,
          START_FEATURE_NUMBER,
          *m_Builder
        );
    }

    //! @brief load a range of data from the dataset, given a particular dataset builder
    //! @param SUBSET the range of data to load
    //! @param FEATURES_STORAGE where to store features that are loaded, must be large enough to hold the subset without resizing
    //! @param RESULTS_STORAGE where to store the corresponding results, must be large enough to hold the subset without resizing
    //! @param START_FEATURE_NUMBER position to store the first feature in FEATURES_STORAGE
    //! @param BUILDER the dataset builder to use
    //! @return # of features actually loaded
    size_t ProteinWithCacheDatasetFromFile::GenerateDataSubsetGivenBuilder
    (
      const math::Range< size_t> &SUBSET,
      linal::MatrixInterface< float> &FEATURES_STORAGE,
      linal::MatrixInterface< float> &RESULTS_STORAGE,
      linal::MatrixInterface< char> &IDS_STORAGE,
      const size_t &START_FEATURE_NUMBER,
      descriptor::DatasetBuilder< biol::AABase> &BUILDER
    )
    {
      pdb::Factory::GetFlagConvertToNaturalAAType()->SetFlag();

      // retrieve the protein ensemble
      util::ShPtrVector< ProteinModelWithCache> models;
      size_t start_feature( 0), total_size( 0);

      BCL_Assert
      (
        BUILDER.GetFeatureSize() == FEATURES_STORAGE.GetNumberCols(),
        "Wrong column size for features " + util::Format()( BUILDER.GetFeatureSize())
        + " storage: " + util::Format()( FEATURES_STORAGE.GetNumberCols())
      );
      BCL_Assert( BUILDER.GetResultSize() == RESULTS_STORAGE.GetNumberCols(), "Wrong column size for results");
      BCL_Assert( BUILDER.GetIdSize() == IDS_STORAGE.GetNumberCols(), "Wrong column size for ids");

      if( BUILDER.GetType().GetDimension() == size_t( 0))
      {
        const math::Range< size_t> closed_range( SUBSET.CloseBorders());
        models = RetrieveEnsemble( closed_range);
        total_size = models.GetSize();
      }
      else
      {
        math::Range< size_t> closed_range( SUBSET.CloseBorders());
        storage::Map< size_t, size_t>::const_iterator itr_lower( m_FeatureStartIdToKey.LowerBound( closed_range.GetMin()));
        storage::Map< size_t, size_t>::const_iterator itr_upper( m_FeatureStartIdToKey.LowerBound( closed_range.GetMax()));
        if( itr_lower == m_FeatureStartIdToKey.End())
        {
          return 0;
        }
        else if( itr_lower->first > closed_range.GetMin())
        {
          --itr_lower;
        }
        if( itr_upper == m_FeatureStartIdToKey.End())
        {
          --itr_upper;
          closed_range = math::Range< size_t>( closed_range.GetMin(), itr_upper->first - 1);
        }
        total_size = closed_range.GetWidth() + 1;
        const math::Range< size_t> model_range
        (
          math::RangeBorders::e_LeftClosed,
          itr_lower->second,
          itr_upper->second,
          math::RangeBorders::e_RightOpen
        );
        models = RetrieveEnsemble( model_range);
        start_feature = closed_range.GetMin() - itr_lower->first;
      }
      iterate::Generic< const descriptor::SequenceInterface< biol::AABase> >
        model_iterator( models.Begin(), models.End());
      descriptor::Dataset dataset( BUILDER( model_iterator, start_feature, total_size));
      // copy the code vectors into the matrices
      size_t data_counter( START_FEATURE_NUMBER);
      std::copy
      (
        dataset.GetFeaturesReference().Begin(),
        dataset.GetFeaturesReference().End(),
        FEATURES_STORAGE[ data_counter]
      );
      // copy the result
      std::copy
      (
        dataset.GetResultsReference().Begin(),
        dataset.GetResultsReference().End(),
        RESULTS_STORAGE[ data_counter]
      );
      // copy the ids
      std::copy
      (
        dataset.GetIdsReference().Begin(),
        dataset.GetIdsReference().End(),
        IDS_STORAGE[ data_counter]
      );

      // return the # of features that were actually loaded
      return dataset.GetSize();
    }

  } // namespace assemble
} // namespace bcl
