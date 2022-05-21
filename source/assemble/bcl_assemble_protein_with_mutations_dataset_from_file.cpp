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
#include "assemble/bcl_assemble_protein_with_mutations_dataset_from_file.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model_with_mutations.h"
#include "biol/bcl_biol_aa_classes.h"
#include "biol/bcl_biol_protein_mutation_set.h"
#include "descriptor/bcl_descriptor_dataset_builder.h"
#include "io/bcl_io_directory_entry.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_serialization.h"
#include "pdb/bcl_pdb_factory.h"
#include "util/bcl_util_wrapper.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {

  //////////
  // data //
  //////////

    const util::SiPtr< const util::ObjectInterface> ProteinWithMutationsDatasetFromFile::s_SequenceInstance
    (
      util::Enumerated< model::RetrieveDataSetBase>::AddInstance( new ProteinWithMutationsDatasetFromFile( false))
    );
    const util::SiPtr< const util::ObjectInterface> ProteinWithMutationsDatasetFromFile::s_ProteinInstance
    (
      util::Enumerated< model::RetrieveDataSetBase>::AddInstance( new ProteinWithMutationsDatasetFromFile( true))
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor, takes the desired AA class and whether or not to require coordinates
    ProteinWithMutationsDatasetFromFile::ProteinWithMutationsDatasetFromFile( const bool &REQUIRE_COORDINATES) :
      m_ProteinStorage( REQUIRE_COORDINATES),
      m_InvertFilter( false),
      m_FractSelfMutations( 0.0),
      m_Builder( new descriptor::DatasetBuilder< biol::Mutation>())
    {
    }

    //! @brief copy constructor, clones builder
    ProteinWithMutationsDatasetFromFile::ProteinWithMutationsDatasetFromFile( const ProteinWithMutationsDatasetFromFile &PARENT) :
      model::RetrieveDataSetBase( PARENT),
      m_ProteinStorage( PARENT.m_ProteinStorage),
      m_FeatureStartIdToKey( PARENT.m_FeatureStartIdToKey),
      m_FilteredMutations( PARENT.m_FilteredMutations),
      m_InvertFilter( PARENT.m_InvertFilter),
      m_FractSelfMutations( PARENT.m_FractSelfMutations),
      m_Builder
      (
        new descriptor::DatasetBuilder< biol::Mutation>
        (
          PARENT.GetFeatureCode(),
          PARENT.GetResultCode(),
          PARENT.GetIdCode()
        )
      )
    {
    }

    //! @brief Clone function
    //! @return pointer to new ProteinWithMutationsDatasetFromFile
    ProteinWithMutationsDatasetFromFile *ProteinWithMutationsDatasetFromFile::Clone() const
    {
      return new ProteinWithMutationsDatasetFromFile( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ProteinWithMutationsDatasetFromFile::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &ProteinWithMutationsDatasetFromFile::GetAlias() const
    {
      static const std::string s_seqname( "SequenceMutationsDirectory"), s_protname( "ProteinMutationsDirectory");
      return m_ProteinStorage.GetRequireCoordinates() ? s_protname : s_seqname;
    }

    //! @brief Set the code / label for the feature (1st part) of the data set
    //! @param CODE the new code
    void ProteinWithMutationsDatasetFromFile::SelectFeatures( const util::ObjectDataLabel &CODE)
    {
      model::RetrieveDataSetBase::SelectFeatures( CODE);
      // this is not strictly necessary, but it allows a user to ask for help over the command line for this retriever
      util::Implementation< descriptor::Base< biol::Mutation, float> > property( GetFeatureCode());
      m_Builder->SetFeatureCode( GetFeatureCode());
    }

    //! @brief Set the code / label for the result (2nd part) of the data set
    //! @param CODE the new code
    void ProteinWithMutationsDatasetFromFile::SelectResults( const util::ObjectDataLabel &CODE)
    {
      model::RetrieveDataSetBase::SelectResults( CODE);
      m_FeatureStartIdToKey.Reset();
      m_Builder->SetResultsCode( GetResultCode());
      const descriptor::Type type( m_Builder->GetType());
      if( type.GetDimension() == size_t( 0))
      {
        // add up the size of each, convoluted by the type
        size_t key_number( 0);
        for( const size_t n_proteins( m_ProteinStorage.GetSize()); key_number < n_proteins; ++key_number)
        {
          m_FeatureStartIdToKey[ key_number] = key_number;
        }
        m_FeatureStartIdToKey[ key_number] = key_number;
      }
      else
      {
        // retrieve the size of each key
        storage::Vector< std::string> all_keys( m_ProteinStorage.GetAllKeys());
        // add up the size of each, convoluted by the type
        size_t total_features( 0), key_number( 0);
        util::GetLogger() << "Computing total # of features\n";
        io::IFStream input;
        for
        (
          storage::Vector< std::string>::const_iterator itr( all_keys.Begin()), itr_end( all_keys.End());
          itr != itr_end;
          ++itr, ++key_number
        )
        {
          io::File::MustOpenIFStream( input, io::File::RemoveLastExtension( *itr) + m_Suffix);
          const size_t n_new_features
          (
            type.GetNumberFeatures( util::StringLineListFromIStream( input).GetSize() * ( m_FractSelfMutations ? 2 : 1))
          );
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
          io::File::CloseClearFStream( input);
        }
        m_FeatureStartIdToKey[ total_features] = key_number;
      }
    }

    //! @brief Set the code / label for the ids (3rd part) of the data set
    //! @param CODE the new code
    void ProteinWithMutationsDatasetFromFile::SelectIds( const util::ObjectDataLabel &CODE)
    {
      model::RetrieveDataSetBase::SelectIds( CODE);
      // this is not strictly necessary, but it allows a user to ask for help over the command line for this retriever
      util::Implementation< descriptor::Base< biol::Mutation, char> > property( GetIdCode());
      m_Builder->SetIdCode( CODE);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief Get the code / label set for the feature (1st part) of the data set with sizes of each property
    //! @return the code / label for the feature (1st part) of the data set with sizes of each property
    //! the feature code set
    model::FeatureLabelSet ProteinWithMutationsDatasetFromFile::GetFeatureLabelsWithSizes() const
    {
      return m_Builder.IsDefined() ? m_Builder->GetFeatureCode().GetLabelsWithSizes() : model::FeatureLabelSet();
    }

    //! @brief Get the code / label for the result (2nd part) of the data set with sizes of each property
    //! @return the code / label for the result (2nd part) of the data set with sizes of each property
    //! the feature code set
    model::FeatureLabelSet ProteinWithMutationsDatasetFromFile::GetResultCodeWithSizes() const
    {
      return m_Builder.IsDefined() ? m_Builder->GetResultCode().GetLabelsWithSizes() : model::FeatureLabelSet();
    }

    //! @brief Get the code / label for the ids of the data set with sizes of each property
    //! @return the code / label for the ids of the data set with sizes of each property
    //! the feature code set
    model::FeatureLabelSet ProteinWithMutationsDatasetFromFile::GetIdCodeWithSizes() const
    {
      return m_Builder.IsDefined() ? m_Builder->GetIdCode().GetLabelsWithSizes() : model::FeatureLabelSet();
    }

    //! @brief generate dataset from a set of ranges
    //! @return generated dataset
    util::ShPtr< descriptor::Dataset> ProteinWithMutationsDatasetFromFile::GenerateDataSet()
    {
      pdb::Factory::GetFlagConvertToNaturalAAType()->SetFlag();

      // determine the total size of the data set
      const size_t total_dataset_size( m_FeatureStartIdToKey.ReverseBegin()->first);

      if( !m_Builder.IsDefined())
      {
        m_Builder =
          util::ShPtr< descriptor::DatasetBuilder< biol::Mutation> >
          (
            new descriptor::DatasetBuilder< biol::Mutation>( GetFeatureCode(), GetResultCode(), GetIdCode())
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

    //! @brief get the number of partitions requested by the user, along with the partition ids
    //! @return the number of partitions requested by the user, along with the partition ids
    storage::Pair< size_t, math::RangeSet< size_t> > ProteinWithMutationsDatasetFromFile::GetNumberPartitionsAndIds() const
    {
      return storage::Pair< size_t, math::RangeSet< size_t> >( 1, math::RangeSet< size_t>( math::Range< size_t>( 0, 0)));
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer ProteinWithMutationsDatasetFromFile::GetSerializer() const
    {
      io::Serializer member_data( m_ProteinStorage.GetSerializer());
      member_data.SetClassDescription
      (
        std::string
        (
          m_ProteinStorage.GetRequireCoordinates()
          ? "Calculates descriptors related to mutations of proteins"
          : "Calculates descriptors related to mutations of protein sequences"
        ) +
        " from a directory or subdirectory"
      );
      member_data.AddInitializer
      (
        "mutation extension",
        "file extension for mutation file. Replacing .pdb/.fasta with this value should give us the mutation file",
        io::Serialization::GetAgent( &m_Suffix),
        ".mutants.txt"
      );
      member_data.AddOptionalInitializer
      (
        "filter",
        "mutations to filter, formatted as protein key to set of mutations to filter, e.g. : (a1q1=(H176G))",
        io::Serialization::GetAgent( &m_FilteredMutations)
      );
      member_data.AddInitializer
      (
        "invert filter",
        "Whether to invert the sense of the filter. If true, only the mutations given in the filter will be returned",
        io::Serialization::GetAgent( &m_InvertFilter),
        "False"
      );

      member_data.AddInitializer
      (
        "add self mutation fraction",
        "Fraction of mutations for which the synonymous mutations should be added to help the model learn that native is native",
        io::Serialization::GetAgentWithRange( &m_FractSelfMutations, 0.0, 1.0),
        "0.0"
      );
      return member_data;
    } // GetParameters

    //! @brief get the nominal (e.g. best estimate without generating the entire dataset) size of the dataset
    //! @return the nominal (e.g. best estimate without generating the entire dataset) size of the dataset
    size_t ProteinWithMutationsDatasetFromFile::GetNominalSize() const
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
    size_t ProteinWithMutationsDatasetFromFile::GenerateDataSubset
    (
      const math::Range< size_t> &SUBSET,
      linal::MatrixInterface< float> &FEATURES_STORAGE,
      linal::MatrixInterface< float> &RESULTS_STORAGE,
      linal::MatrixInterface< char> &IDS_STORAGE,
      const size_t &START_FEATURE_NUMBER
    )
    {
      if( !m_Builder.IsDefined())
      {
        m_Builder =
          util::ShPtr< descriptor::DatasetBuilder< biol::Mutation> >
          (
            new descriptor::DatasetBuilder< biol::Mutation>( GetFeatureCode(), GetResultCode(), GetIdCode())
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
    size_t ProteinWithMutationsDatasetFromFile::GenerateDataSubsetGivenBuilder
    (
      const math::Range< size_t> &SUBSET,
      linal::MatrixInterface< float> &FEATURES_STORAGE,
      linal::MatrixInterface< float> &RESULTS_STORAGE,
      linal::MatrixInterface< char> &IDS_STORAGE,
      const size_t &START_FEATURE_NUMBER,
      descriptor::DatasetBuilder< biol::Mutation> &BUILDER
    )
    {
      // retrieve the protein ensemble
      util::ShPtrVector< biol::ProteinMutationSet> models;
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
        auto temp_models( m_ProteinStorage.RetrieveEnsemble( closed_range));
        models.AllocateMemory( temp_models.GetSize());
        for
        (
          auto itr_tmp_models( temp_models.Begin()), itr_tmp_models_end( temp_models.End());
          itr_tmp_models != itr_tmp_models_end;
          ++itr_tmp_models
        )
        {
          models.PushBack
          (
            util::ShPtr< biol::ProteinMutationSet>
            (
              new biol::ProteinMutationSet
              (
                **itr_tmp_models,
                m_ProteinStorage.GetRequireCoordinates(),
                GetMutationsFromProteinModel( **itr_tmp_models)
              )
            )
          );
        }
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
        auto temp_models( m_ProteinStorage.RetrieveEnsemble( model_range));
        models.AllocateMemory( temp_models.GetSize());
        for
        (
          auto itr_tmp_models( temp_models.Begin()), itr_tmp_models_end( temp_models.End());
          itr_tmp_models != itr_tmp_models_end;
          ++itr_tmp_models
        )
        {
          models.PushBack
          (
            util::ShPtr< biol::ProteinMutationSet>
            (
              new biol::ProteinMutationSet
              (
                **itr_tmp_models,
                m_ProteinStorage.GetRequireCoordinates(),
                GetMutationsFromProteinModel( **itr_tmp_models)
              )
            )
          );
        }

        start_feature = closed_range.GetMin() - itr_lower->first;
      }
      iterate::Generic< const descriptor::SequenceInterface< biol::Mutation> >
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

    storage::Vector< biol::Mutation> ProteinWithMutationsDatasetFromFile::GetMutationsFromProteinModel
    (
      const ProteinModel &MODEL
    ) const
    {
      // Get the filename
      util::ShPtr< util::Wrapper< std::string> > sp_filename_wrapper
      (
        MODEL.GetProteinModelData()->GetData( ProteinModelData::e_PDBFile)
      );

      // Remove the last extension
      std::string basename
      (
        io::File::RemoveLastExtension( io::File::RemoveCompressionExtension( sp_filename_wrapper->GetData()))
      );
      io::DirectoryEntry entry( basename + m_Suffix);
      storage::Vector< biol::Mutation> mutations;
      io::IFStream input;
      io::File::MustOpenIFStream( input, basename + m_Suffix);
      if( entry.DoesExist())
      {
        storage::Set< biol::Mutation> empty_set;
        auto itr_filter_map( m_FilteredMutations.Find( io::File::RemovePath( basename)));
        const storage::Set< biol::Mutation> &filtered_mutations
        (
          itr_filter_map == m_FilteredMutations.End()
          ? empty_set
          : itr_filter_map->second
        );
        storage::Vector< std::string> mutation_strs( util::StringListFromIStream( input));
        for( auto itr( mutation_strs.Begin()), itr_end( mutation_strs.End()); itr != itr_end; ++itr)
        {
          biol::Mutation candidate_mutation( biol::Mutation::FromString( util::SplitString( util::TrimString( *itr))( 0)));
          if( filtered_mutations.Contains( candidate_mutation) == m_InvertFilter)
          {
            mutations.PushBack( candidate_mutation);
            if( m_FractSelfMutations && random::GetGlobalRandom().Double() < m_FractSelfMutations)
            {
              biol::Mutation self_mut
              (
                candidate_mutation.GetResidueNumber(),
                candidate_mutation.GetNativeType(),
                candidate_mutation.GetNativeType()
              );
              if( mutations.Find( self_mut) >= mutations.GetSize())
              {
                mutations.PushBack( self_mut);
              }
            }
          }
        }
      }
      else
      {
        BCL_MessageCrt( "Mutations file: " + basename + m_Suffix + " not found, skipping!");
      }
      io::File::CloseClearFStream( input);
      return mutations;
    }

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERR_STREAM stream to write out errors to
    bool ProteinWithMutationsDatasetFromFile::ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERR_STREAM)
    {
      return m_ProteinStorage.ReadInitializerSuccessHook( LABEL, ERR_STREAM);
    }

  } // namespace assemble
} // namespace bcl
