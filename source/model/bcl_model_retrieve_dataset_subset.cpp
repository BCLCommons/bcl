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
#include "model/bcl_model_retrieve_dataset_subset.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_binary_serialize.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_serialization.h"
#include "sched/bcl_sched_scheduler_interface.h"
#include "sched/bcl_sched_thunk_job.h"

// external includes - sorted alphabetically
#include "bcl_version.h"
#include <fstream>

namespace bcl
{
  namespace model
  {

  //////////
  // data //
  //////////

    // add the interface to the set of known implementations
    const util::SiPtr< const util::ObjectInterface> RetrieveDatasetSubset::s_Instance
    (
      util::Enumerated< RetrieveDataSetBase>::AddInstance( new RetrieveDatasetSubset())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    //! @param FILENAME the file to retrieve the data set from
    //! @param NUMBER_CHUNKS # of chunks the file should be split into (conceptually), used to divide the file into disparate datasets
    //! @param CHUNKS chunks to load
    RetrieveDatasetSubset::RetrieveDatasetSubset
    (
      const std::string &FILENAME,
      const size_t &NUMBER_CHUNKS,
      const math::RangeSet< size_t> &CHUNKS
    ) :
      m_Filename( FILENAME),
      m_NumberChunks( NUMBER_CHUNKS),
      m_ChunkRanges( CHUNKS)
    {
    }

    //! @brief Clone is the virtual Clone constructor
    //! @return a pointer to new LocatorAtomFurthest which is a copy of this
    RetrieveDatasetSubset *RetrieveDatasetSubset::Clone() const
    {
      return new RetrieveDatasetSubset( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &RetrieveDatasetSubset::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &RetrieveDatasetSubset::GetAlias() const
    {
      static const std::string s_Name( "Subset");
      return s_Name;
    }

    //! @brief Get the code / label set for the feature (1st part) of the data set with sizes of each property
    //! @return the code / label for the feature (1st part) of the data set with sizes of each property
    //! the feature code set
    FeatureLabelSet RetrieveDatasetSubset::GetFeatureLabelsWithSizes() const
    {
      std::ifstream input;
      input.open( m_Filename.c_str(), std::ios::binary | std::ios::in);
      FeatureLabelSet all_feature_labels( ReadFeatureLabelSetBinary( input));
      io::File::CloseClearFStream( input);

      // determine which features were actually selected
      storage::Vector< size_t> feature_sizes;

      // check; if no features were given, then set all features
      if( GetFeatureCode().IsScalar())
      {
        return all_feature_labels;
      }
      feature_sizes.AllocateMemory( GetFeatureCode().GetNumberArguments());

      // determine the size of each result that was selected
      for
      (
        util::ObjectDataLabel::const_iterator
          itr_label( GetFeatureCode().Begin()),
          itr_label_end( GetFeatureCode().End());
        itr_label != itr_label_end;
        ++itr_label
      )
      {
        // find this feature label in the feature_labels object
        feature_sizes.PushBack( all_feature_labels.GetPropertyIndices( *itr_label).GetSize());
      }

      return FeatureLabelSet( all_feature_labels.GetAlias(), GetFeatureCode().GetArguments(), feature_sizes);
    }

    //! @brief Get the code / label for the result (2nd part) of the data set with sizes of each property
    //! @return the code / label for the result (2nd part) of the data set with sizes of each property
    //! the feature code set
    FeatureLabelSet RetrieveDatasetSubset::GetResultCodeWithSizes() const
    {
      std::ifstream input;
      input.open( m_Filename.c_str(), std::ios::binary | std::ios::in);
      // skip the feature labels
      ReadFeatureLabelSetBinary( input);
      // read the result labels
      FeatureLabelSet all_result_labels( ReadFeatureLabelSetBinary( input));
      io::File::CloseClearFStream( input);

      // check; if no result code was given, then return the entire result code object
      if( GetResultCode().IsScalar())
      {
        return all_result_labels;
      }

      // determine which results were actually selected
      storage::Vector< size_t> results_sizes;
      results_sizes.AllocateMemory( GetResultCode().GetNumberArguments());

      // determine the size of each result that was selected
      for
      (
        util::ObjectDataLabel::const_iterator
          itr_label( GetResultCode().Begin()),
          itr_label_end( GetResultCode().End());
        itr_label != itr_label_end;
        ++itr_label
      )
      {
        // find this feature label in the feature_labels object
        results_sizes.PushBack( all_result_labels.GetPropertyIndices( *itr_label).GetSize());
      }

      return FeatureLabelSet( all_result_labels.GetAlias(), GetResultCode().GetArguments(), results_sizes);
    }

    //! @brief Get the code / label for the ids of the data set with sizes of each property
    //! @return the code / label for the ids of the data set with sizes of each property
    //! the feature code set
    FeatureLabelSet RetrieveDatasetSubset::GetIdCodeWithSizes() const
    {
      std::ifstream input;
      input.open( m_Filename.c_str(), std::ios::binary | std::ios::in);
      // skip the feature labels
      ReadFeatureLabelSetBinary( input);
      // skip the result labels
      ReadFeatureLabelSetBinary( input);

      // read the id labels
      FeatureLabelSet all_id_labels( ReadFeatureLabelSetBinary( input));
      io::File::CloseClearFStream( input);

      // check; if no id code was given, then return the entire result code object
      if( GetIdCode().IsScalar())
      {
        return all_id_labels;
      }

      // determine which results were actually selected
      storage::Vector< size_t> id_sizes;
      id_sizes.AllocateMemory( GetIdCode().GetNumberArguments());

      // determine the size of each result that was selected
      for
      (
        util::ObjectDataLabel::const_iterator itr_label( GetIdCode().Begin()), itr_label_end( GetIdCode().End());
        itr_label != itr_label_end;
        ++itr_label
      )
      {
        // find this feature label in the feature_labels object
        id_sizes.PushBack( all_id_labels.GetPropertyIndices( *itr_label).GetSize());
      }

      return FeatureLabelSet( all_id_labels.GetAlias(), GetIdCode().GetArguments(), id_sizes);
    }

    //! @brief get the number of partitions requested by the user, along with the partition ids
    //! @return the number of partitions requested by the user, along with the partition ids
    storage::Pair< size_t, math::RangeSet< size_t> > RetrieveDatasetSubset::GetNumberPartitionsAndIds() const
    {
      return storage::Pair< size_t, math::RangeSet< size_t> >( m_NumberChunks, m_ChunkRanges);
    }

    //! @brief write attributes for the arff file header
    //! @param OSTREAM output stream
    //! @param ID_LABELS labels for IDS
    //! @param FEATURE_LABELS labels for the features
    //! @param RESULT_LABELS labels for the results
    void WriteArffAttributes
    (
      std::ostream &OSTREAM,
      const FeatureLabelSet &ID_LABELS,
      const FeatureLabelSet &FEATURE_LABELS,
      const FeatureLabelSet &RESULT_LABELS
    )
    {
      storage::Vector< std::string> ids( ID_LABELS.GetMemberLabels().Begin(), ID_LABELS.GetMemberLabels().End());

      // get the labels; one per actual numeric column
      FeatureLabelSet features_split( FEATURE_LABELS.SplitFeatureLabelSet( true));
      FeatureLabelSet results_split( RESULT_LABELS.SplitFeatureLabelSet( true));
      storage::Vector< std::string> features
      (
        features_split.GetMemberLabels().Begin(),
        features_split.GetMemberLabels().End()
      );
      storage::Vector< std::string> results
      (
        results_split.GetMemberLabels().Begin(),
        results_split.GetMemberLabels().End()
      );
      util::SiPtrVector< storage::Vector< std::string> > svect
      (
        util::SiPtrVector< storage::Vector< std::string> >::Create
        (
          ids,
          features,
          results
        )
      );
      // change all double quotes to single for weka
      for
      (
        util::SiPtrVector< storage::Vector< std::string> >::iterator itr_vec( svect.Begin()), itr_vec_end( svect.End());
        itr_vec != itr_vec_end;
        ++itr_vec
      )
      {
        for
        (
          storage::Vector< std::string>::iterator itr( ( *itr_vec)->Begin()), itr_end( ( *itr_vec)->End());
          itr != itr_end;
          ++itr
        )
        {
          std::replace( itr->begin(), itr->end(), '"', '\'');
        }
      }
      if( ids.GetSize())
      {
        OSTREAM << "% IDs section (used only for interpretation)\n";
        OSTREAM << "@attribute \"" << util::Join( "\" string\n@attribute \"", ids) << "\" string\n";
      }
      if( features.GetSize())
      {
        OSTREAM << "% Features section (descriptors used to predict with)\n";
        OSTREAM << "@attribute \"" << util::Join( "\" numeric\n@attribute \"", features) << "\" numeric\n";
      }
      if( results.GetSize())
      {
        OSTREAM << "% Results. For classification problems, replace \"numeric\" type below with the allowed class values\n";
        OSTREAM << "@attribute \"" << util::Join( "\" numeric\n@attribute \"", results) << "\" numeric\n";
      }
      OSTREAM << "\n@data\n";
    }

    //! @brief helper function to write out id labels with commas between the individual properties
    //! @param OSTREAM output stream
    //! @param ID_ROW row of the id matrix
    //! @param SIZES sizes of each descriptor in the row
    void WriteIdLabelsCSV( std::ostream &OSTREAM, const char *const &ID_ROW, const storage::Vector< size_t> &SIZES)
    {
      if( !SIZES.IsEmpty())
      {
        // write the first property to the ostream
        OSTREAM.write( ID_ROW, SIZES.FirstElement());

        // start of the next descriptor, updated in for loop
        const char *descriptor_start( ID_ROW + SIZES.FirstElement());

        for
        (
          storage::Vector< size_t>::const_iterator itr( SIZES.Begin() + 1), itr_end( SIZES.End());
          itr != itr_end;
          ++itr
        )
        {
          // write a comma
          OSTREAM << ',';

          // write the given property to the ostream
          OSTREAM.write( descriptor_start, *itr);

          // update the start pointer for the next descriptor
          descriptor_start += *itr;
        }
      }
    }

    //! @brief helper function to write features, results, and ids to a csv file
    void WriteCSV
    (
      std::ostream &OUTPUT,
      const linal::Matrix< float> &FEATURES,
      const linal::Matrix< float> &RESULTS,
      const linal::Matrix< char> &IDS,
      const storage::Vector< size_t> &ID_SIZES,
      const size_t &NUMBER_ROWS,
      const bool &ALLOW_INCOMPLETE_RECORDS
    )
    {
      // determine the feature and result sizes
      const size_t feature_size( FEATURES.GetNumberCols()), result_size( RESULTS.GetNumberCols());
      const size_t id_size( IDS.GetNumberCols());

      // csv output
      // write out each feature/result on one line, with line size exactly equal to line_width
      for( size_t feature_number( 0); feature_number < NUMBER_ROWS; ++feature_number)
      {
        if( !ALLOW_INCOMPLETE_RECORDS && !RESULTS.GetRow( feature_number).IsDefined())
        {
          continue;
        }

        // write out the ids
        WriteIdLabelsCSV( OUTPUT, IDS[ feature_number], ID_SIZES);

        // keep track of whether anything has been written on this line; necessary to keep all commas internal
        bool have_written_first( id_size);
        const float *itr_feature( FEATURES[ feature_number]), *itr_feature_end( itr_feature + feature_size);
        if( !have_written_first && feature_size)
        {
          have_written_first = true;
          // write out the first value
          io::Serialize::Write( *itr_feature, OUTPUT);
          ++itr_feature;
        }
        // write out the features
        for( ; itr_feature != itr_feature_end; ++itr_feature)
        {
          // csv
          OUTPUT << ',';
          io::Serialize::Write( *itr_feature, OUTPUT);
        }

        const float *itr_result( RESULTS[ feature_number]), *itr_result_end( itr_result + result_size);
        if( !have_written_first && result_size)
        {
          // write out the first value
          io::Serialize::Write( *itr_result, OUTPUT);
          ++itr_result;
        }
        // write out the results
        for( ; itr_result != itr_result_end; ++itr_result)
        {
          OUTPUT << ',';
          io::Serialize::Write( *itr_result, OUTPUT);
        }
        OUTPUT << '\n';
      }
    }

    //! @brief helper function to write features, results, and ids to a csv file
    void WriteBin
    (
      std::ostream &OUTPUT,
      const linal::Matrix< float> &FEATURES,
      const linal::Matrix< float> &RESULTS,
      const linal::Matrix< char> &IDS,
      const size_t &NUMBER_ROWS,
      const bool &ALLOW_INCOMPLETE_RECORDS
    )
    {
      // determine the feature and result sizes
      const size_t feature_size( FEATURES.GetNumberCols()), result_size( RESULTS.GetNumberCols());
      const size_t id_size( IDS.GetNumberCols());

      // write out each feature/result on one line, with line size exactly equal to line_width
      for( size_t feature_number( 0); feature_number < NUMBER_ROWS; ++feature_number)
      {
        if( !ALLOW_INCOMPLETE_RECORDS && !RESULTS.GetRow( feature_number).IsDefined())
        {
          continue;
        }
        // write out the ids
        for
        (
          const char *itr_feature( IDS[ feature_number]), *itr_feature_end( itr_feature + id_size);
          itr_feature != itr_feature_end;
          ++itr_feature
        )
        {
          // write binary output
          io::BinarySerialize::Write( *itr_feature, OUTPUT);
        }

        // write out the features
        for
        (
          const float *itr_feature( FEATURES[ feature_number]), *itr_feature_end( itr_feature + feature_size);
          itr_feature != itr_feature_end;
          ++itr_feature
        )
        {
          // binary
          io::BinarySerialize::Write( *itr_feature, OUTPUT);
        }

        // write out the results
        for
        (
          const float *itr_feature( RESULTS[ feature_number]), *itr_feature_end( itr_feature + result_size);
          itr_feature != itr_feature_end;
          ++itr_feature
        )
        {
          // write binary output
          io::BinarySerialize::Write( *itr_feature, OUTPUT);
        }
      }
    }

    // local enum to specify a storage type
    enum StorageType
    {
      e_Bin,  // .bin file
      e_Csv,  // .csv file, or csv.bz2, etc.
      e_Arff  // .arff file, used for interaction with WEKA (http://www.cs.waikato.ac.nz/ml/weka/arff.html)
    };

    //! @class struct to allow for threaded dataset creation
    struct BCL_API StoreThreadData
    {
      math::Range< size_t> m_Range; //!< Dataset range to construct
      util::Implementation< RetrieveDataSetBase> m_Retriever; //!< Retriever of the dataset
      linal::Matrix< float> m_Features; //!< Features calculated by the thread
      linal::Matrix< float> m_Results;  //!< Results calculated by the thread
      linal::Matrix< char>  m_Ids;      //!< IDs calculated by the thread
      StorageType           m_StorageType; //!< Type of storage to use
      size_t                m_NumberFeaturesComputed; //!< Actual # of features computed
      storage::Vector< size_t> m_IdSizes;
      std::string           m_LastString;
      bool                  m_AllowIncompleteRecords;

      //! @brief default constructor
      StoreThreadData() : m_NumberFeaturesComputed( 0) {}

      //! function that actually works with the thread
      void RunThread()
      {
        std::stringstream stream( std::ios::out | ( m_StorageType == e_Bin ? std::ios::binary : std::ios::out));
        m_NumberFeaturesComputed = m_Retriever->GenerateDataSubset( m_Range, m_Features, m_Results, m_Ids, 0);
        if( m_StorageType == e_Bin)
        {
          WriteBin( stream, m_Features, m_Results, m_Ids, m_NumberFeaturesComputed, m_AllowIncompleteRecords);
        }
        else
        {
          WriteCSV( stream, m_Features, m_Results, m_Ids, m_IdSizes, m_NumberFeaturesComputed, m_AllowIncompleteRecords);
        }
        m_LastString = stream.str();
      }
    };

    // trivial helper functions to allow usage of StoreThreadData in the scheduling framework
    std::ostream &operator <<( std::ostream &STREAM, const StoreThreadData &)
    {
      return STREAM;
    }

    // trivial helper functions to allow usage of StoreThreadData in the scheduling framework
    std::istream &operator >>( std::istream &STREAM, StoreThreadData &)
    {
      return STREAM;
    }

    //! @brief store the data from a dataset into a master dataset
    //! @param FILENAME file to store the dataset in
    //! @param DATA_SET_RETRIEVER object used to retrieve the dataset
    //! @param BLOCK_SIZE_MB maximum MB of data to generate per thread before writing to disk
    //! @param ALLOW_INCOMPLETE_RECORDS true to allow incomplete results to be stored; by default, all results
    //!        must be given for each record
    void RetrieveDatasetSubset::StoreMasterDataset
    (
      const std::string &FILENAME,
      RetrieveDataSetBase &DATA_SET_RETRIEVER,
      const double &BLOCK_SIZE_MB,
      const bool &ALLOW_INCOMPLETE_RECORDS
    )
    {
      // init storage type, assume binary file, since that is the default format
      StorageType storage_type( e_Bin);

      // look for csv in the filename
      const size_t csv_pos( FILENAME.rfind( ".csv"));
      // if there is a csv extension, assume .csv format, but still allow compressed alternatives
      if( csv_pos != std::string::npos && ( csv_pos + 4 == FILENAME.size() || FILENAME[ csv_pos + 4] == '.'))
      {
        storage_type = e_Csv;
      }
      else
      {
        // check for arff format
        const size_t arff_pos( FILENAME.rfind( ".arff"));
        if( arff_pos != std::string::npos && ( arff_pos + 5 == FILENAME.size() || FILENAME[ arff_pos + 5] == '.'))
        {
          storage_type = e_Arff;
        }
      }
      io::OFStream output;

      if( storage_type == e_Bin)
      {
        io::File::MustOpenOFStream( output, FILENAME, std::ios::binary);
      }
      else
      {
        io::File::MustOpenOFStream( output, FILENAME);
      }
      BCL_Assert( output.is_open(), "Could not open " + FILENAME + " for writing");

      // get the feature and result labels with sizes
      FeatureLabelSet feature_labels( DATA_SET_RETRIEVER.GetFeatureLabelsWithSizes());
      FeatureLabelSet result_labels( DATA_SET_RETRIEVER.GetResultCodeWithSizes());
      FeatureLabelSet id_labels( DATA_SET_RETRIEVER.GetIdCodeWithSizes());

      if( storage_type == e_Bin)
      {
        // write out the feature and result labels
        WriteFeatureLabelSetBinary( output, feature_labels);
        WriteFeatureLabelSetBinary( output, result_labels);
        WriteFeatureLabelSetBinary( output, id_labels);
      }
      else if( storage_type == e_Arff)
      {
        // write out the arff headers
        output << "% Written by " << GetVersion().GetDetailedInfoString() << '\n';
        std::string dataset_retriever( DATA_SET_RETRIEVER.GetCompleteSerializer().GetLabel().ToString());
        std::replace( dataset_retriever.begin(), dataset_retriever.end(), '"', '\'');
        output << "@relation \"" << dataset_retriever << "\"\n\n";
        WriteArffAttributes
        (
          output,
          id_labels,
          feature_labels,
          result_labels
        );
      }

      // determine the nominal size of the entire dataset
      const size_t dataset_size( DATA_SET_RETRIEVER.GetNominalSize());

      BCL_MessageStd
      (
        "Values per feature: " + util::Format()( feature_labels.GetSize())
        + " Values per result: " + util::Format()( result_labels.GetSize())
        + " Characters per id: " + util::Format()( id_labels.GetSize())
        + " Nominal dataset # feature rows: " + util::Format()( dataset_size)
      );

      // determine the feature and result sizes
      const size_t feature_size( feature_labels.GetSize()), result_size( result_labels.GetSize());
      const size_t id_size( id_labels.GetSize());

      // determine the # of features and results put together
      const size_t number_features_and_results( sizeof( float) * ( feature_size + result_size) + id_size);

      size_t number_features_to_generate_at_a_time( 0);

      // determine # cpus
      const size_t n_cpus( sched::GetScheduler().GetNumberUnusedCPUS());

      // determine the number of features to generate at a time
      if( DATA_SET_RETRIEVER.SupportsEfficientSubsetLoading())
      {
        const size_t max_desired_values_at_a_time( ( 1 << 20) * BLOCK_SIZE_MB); // generate <= 32 MB worth of floats at a time
        number_features_to_generate_at_a_time =
          std::min( 1 + max_desired_values_at_a_time / number_features_and_results, dataset_size);
        if( n_cpus > size_t( 1) && dataset_size < max_desired_values_at_a_time / number_features_and_results * n_cpus)
        {
          number_features_to_generate_at_a_time = size_t( 1) + dataset_size / n_cpus;
        }
      }
      else
      {
        BCL_MessageCrt
        (
          "Your choice of dataset retrievers requires that the entire dataset be generated at once\n"
          + std::string( "if this exceeds memory resources, use a different means of retrieval")
        );
        // have to load the entire dataset
        number_features_to_generate_at_a_time = dataset_size;
      }

      // temp_features, temp_results will store features and results from each round, until they can be written out
      // to avoid reallocating memory, just create it once
      linal::Matrix< float> temp_features( number_features_to_generate_at_a_time, feature_size);
      linal::Matrix< float> temp_results( number_features_to_generate_at_a_time, result_size);
      linal::Matrix< char> temp_ids( number_features_to_generate_at_a_time, id_size);
      const storage::Vector< size_t> &id_sizes( id_labels.GetPropertySizes());

      // store the final dataset size
      size_t final_dataset_size( 0);

      // handle the serial case
      if
      (
        number_features_to_generate_at_a_time == dataset_size
        || sched::GetScheduler().GetNumberUnusedCPUS() <= size_t( 1)
      )
      {
        // generate the features in several rounds, if possible
        // to avoid the memory overhead of holding the entire dataset at once in memory
        for
        (
          size_t lower_limit( 0);
          lower_limit < dataset_size;
          lower_limit += number_features_to_generate_at_a_time
        )
        {
          // get the next range to extract from the dataset
          math::Range< size_t> next_dataset_range
          (
            math::RangeBorders::e_LeftClosed,
            lower_limit,
            lower_limit + number_features_to_generate_at_a_time,
            math::RangeBorders::e_RightOpen
          );

          // generate the features; store how many were actually generated
          const size_t actual_number_of_generated_features
          (
            DATA_SET_RETRIEVER.GenerateDataSubset
            (
              next_dataset_range,
              temp_features,
              temp_results,
              temp_ids,
              0
            )
          );

          final_dataset_size += actual_number_of_generated_features;

          std::stringstream output_str;
          if( storage_type == e_Bin)
          {
            WriteBin( output_str, temp_features, temp_results, temp_ids, actual_number_of_generated_features, ALLOW_INCOMPLETE_RECORDS);
          }
          else
          {
            WriteCSV( output_str, temp_features, temp_results, temp_ids, id_sizes, actual_number_of_generated_features, ALLOW_INCOMPLETE_RECORDS);
          }
          const std::string output_str_out( output_str.str());
          output.write( output_str_out.data(), output_str_out.size());
        }
      }
      else
      {
        BCL_MessageStd
        (
          "Computing " + util::Format()( number_features_to_generate_at_a_time) + " features per thread launch"
        );
        // create ranges for each thread
        const size_t expected_thread_launches
        (
          ( dataset_size + number_features_to_generate_at_a_time - 1)
          / number_features_to_generate_at_a_time
        );
        const size_t n_threads( std::max( size_t( 1), std::min( expected_thread_launches, sched::GetNumberCPUs())));
        std::vector< StoreThreadData> thread_data( n_threads);
        util::ShPtrVector< sched::JobInterface> jobs( n_threads);

        // initialize data for each thread
        for( size_t thread_number( 0); thread_number < n_threads; ++thread_number)
        {
          thread_data[ thread_number].m_Features = temp_features;
          thread_data[ thread_number].m_Results = temp_results;
          thread_data[ thread_number].m_Ids = temp_ids;
          thread_data[ thread_number].m_Retriever = util::Implementation< RetrieveDataSetBase>( DATA_SET_RETRIEVER);
          thread_data[ thread_number].m_IdSizes = id_sizes;
          thread_data[ thread_number].m_StorageType = storage_type;
          thread_data[ thread_number].m_AllowIncompleteRecords = ALLOW_INCOMPLETE_RECORDS;
          jobs( thread_number) =
            util::ShPtr< sched::JobInterface>
            (
              new sched::ThunkJob< StoreThreadData, void>( 0, thread_data[ thread_number], &StoreThreadData::RunThread)
            );
        }
        size_t threads_last_round( 0);
        size_t threads_this_round( 0);
        size_t lower_limit( 0);
        for
        (
          ;
          threads_this_round < n_threads && lower_limit < dataset_size;
          ++threads_this_round, lower_limit += number_features_to_generate_at_a_time
        )
        {
          thread_data[ threads_this_round].m_Range =
            math::Range< size_t>
            (
              math::RangeBorders::e_LeftClosed,
              lower_limit,
              lower_limit + number_features_to_generate_at_a_time,
              math::RangeBorders::e_RightOpen
            );
          sched::GetScheduler().SubmitJob( jobs( threads_this_round));
        }

        size_t n_joined( 0);

        // status bar length
        const size_t status_bar_length( 20);

        while( lower_limit < dataset_size || threads_this_round)
        {
          threads_last_round = threads_this_round;
          threads_this_round = 0;
          for( size_t thread_number( 0); thread_number < threads_last_round; ++thread_number)
          {
            sched::GetScheduler().Join( jobs( thread_number));
            ++n_joined;

            // get a reference to the thread data
            const StoreThreadData &thread_data_reference( thread_data[ thread_number]);
            const size_t actual_number_of_generated_features( thread_data_reference.m_NumberFeaturesComputed);
            final_dataset_size += actual_number_of_generated_features;

            // copy the thread's data so that it can immediately be lauched on the next block while this block is
            const std::string thread_chars( thread_data_reference.m_LastString);

            // assign the thread a new work block
            if( lower_limit < dataset_size)
            {
              thread_data[ threads_this_round].m_Range =
                math::Range< size_t>
                (
                  math::RangeBorders::e_LeftClosed,
                  lower_limit,
                  lower_limit + number_features_to_generate_at_a_time,
                  math::RangeBorders::e_RightOpen
                );
              sched::GetScheduler().SubmitJob( jobs( threads_this_round));
              lower_limit += number_features_to_generate_at_a_time;
              ++threads_this_round;
            }

            // determine approximate # of features generated
            const size_t n_features_gen( std::min( n_joined * number_features_to_generate_at_a_time, dataset_size));

            // determine progress percent
            const size_t percent( float( n_joined) * 100.0f / float( expected_thread_launches));

            // determine number of stars in the status bar
            const size_t number_stars( percent * status_bar_length / 100);

            const std::string status
            (
              "["
              + std::string( number_stars, '*')
              + std::string( status_bar_length - number_stars, ' ')
              + "] "
              + util::Format()( percent) + "% "
              + util::Format()( final_dataset_size) + " / "
              + util::Format()( dataset_size - n_features_gen + final_dataset_size)
              + " features generated"
            );
            util::GetLogger().LogStatus( status);

            // write out the data
            output.write( thread_chars.data(), thread_chars.size());
          }
        }
      }

      BCL_MessageStd( "Actual final dataset size: " + util::Format()( final_dataset_size))

      // close the file
      output.close();
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief generate dataset
    //! @return generated dataset
    util::ShPtr< descriptor::Dataset>
    RetrieveDatasetSubset::GenerateDataSet()
    {
      if( m_NumberChunks > 1)
      {
        // ensure the range set does not exceed the # of ranges
        BCL_Assert
        (
          m_ChunkRanges.GetMax() < m_NumberChunks,
          GetClassIdentifier() + " received range set with max >= # ranges! # ranges = " + util::Format()( m_NumberChunks)
          + " " + m_ChunkRanges.AsString()
        );

        // ensure the range set is not empty
        BCL_Assert
        (
          !m_ChunkRanges.IsEmpty(),
          GetClassIdentifier() + " received empty range set"
        );
      }
      else
      {
        m_NumberChunks = 1;
        m_ChunkRanges = math::RangeSet< size_t>( math::Range< size_t>( 0, 0));
      }

      // BCL_MessageStd( "total/actual size: " + util::Format()( total_size) + "/" + util::Format()( actual_dataset_size));

      // read the data set from the file
      std::ifstream input;

      // open the file, determine the size, setup the feature and result columns
      const size_t total_size( GetTotalSizeAndParseLabels( input));

      // get the position before reading this line
      std::istream::pos_type first_line_position( input.tellg());

      // determine the total # of features/results on each line
      const size_t total_features_results( m_Features.GetInputFeatureSize());

      // compute bytes per feautre
      std::istream::pos_type bytes_per_feature( sizeof( float) * total_features_results + m_Ids.GetInputFeatureSize());

      // get the actual ranges that will be taken from the dataset
      storage::Vector< math::Range< size_t> > ranges;
      ranges.AllocateMemory( m_ChunkRanges.GetRanges().GetSize());

      size_t actual_dataset_size( 0); // track the actual size of the resulting dataset
      for
      (
        storage::Set< math::Range< size_t> >::const_iterator
          itr_range( m_ChunkRanges.GetRanges().Begin()), itr_range_end( m_ChunkRanges.GetRanges().End());
        itr_range != itr_range_end;
        ++itr_range
      )
      {
        ranges.PushBack
        (
          RetrieveDataSetBase::GetStartEndPositionOfRange( *itr_range, m_NumberChunks, total_size)
        );
        actual_dataset_size += ranges.LastElement().GetWidth();
      }

      BCL_MessageStd
      (
        "result columns: " + util::Format()( m_Results.GetColumnIndices()( 0))
        + " - " + util::Format()( m_Results.GetColumnIndices().LastElement())
      );

      // create a new dataset
      util::ShPtr< descriptor::Dataset> dataset
      (
        new descriptor::Dataset
        (
          actual_dataset_size,
          GetFeatureLabelsWithSizes(),
          GetResultCodeWithSizes(),
          GetIdCodeWithSizes()
        )
      );

      // get referencest to the underlying matrices
      linal::MatrixReference< float> loaded_features( dataset->GetFeaturesReference());
      linal::MatrixReference< float> loaded_results( dataset->GetResultsReference());
      linal::MatrixReference< char> loaded_ids( dataset->GetIdsReference());

      // create a linal::Vector< float> object to hold the data loaded in
      linal::Vector< float> temp_vector( total_features_results);
      linal::Vector< char> temp_id_vector( loaded_ids.GetNumberCols());
      const size_t id_size( temp_id_vector.GetSize());

      // keep track of how many rows have been loaded
      size_t loaded_feature_counter( 0);

      // walk over the lines desired, loading only those columns and results that we want
      for
      (
        storage::Vector< math::Range< size_t> >::const_iterator itr_range( ranges.Begin()), itr_range_end( ranges.End());
        itr_range != itr_range_end;
        ++itr_range
      )
      {
        // move to the appropriate line
        input.seekg( first_line_position + std::istream::pos_type( itr_range->GetMin()) * bytes_per_feature);

        BCL_MessageVrb( "stream position: " + util::Format()( input.tellg()));

        // load in the desired lines
        for
        (
          size_t line_number( 0), number_lines( itr_range->GetWidth());
          line_number < number_lines;
          ++line_number, ++loaded_feature_counter
        )
        {
          // load the ids into temp_vector
          for( size_t feature_number( 0); feature_number < id_size; ++feature_number)
          {
            io::BinarySerialize::Read( temp_id_vector( feature_number), input);
          }

          // load the features and results into temp_vector
          for( size_t feature_number( 0); feature_number < total_features_results; ++feature_number)
          {
            io::BinarySerialize::Read( temp_vector( feature_number), input);
          }

          // put the chosen columns into the appropriate matrices
          m_Ids( temp_id_vector, loaded_ids[ loaded_feature_counter]);
          m_Features( temp_vector, loaded_features[ loaded_feature_counter]);
          m_Results( temp_vector, loaded_results[ loaded_feature_counter]);
        }
      }

      io::File::CloseClearFStream( input);

      // return the generated data set
      return dataset;
    } // GenerateDataSet

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer RetrieveDatasetSubset::GetSerializer() const
    {
      io::Serializer member_data;
      member_data.SetClassDescription( "reads a .bin file created with GenerateDataset");

      member_data.AddInitializer
      (
        "filename",
        "name of the file from which to load the dataset",
        io::Serialization::GetAgentInputFilename( &m_Filename)
      );

      member_data.AddInitializer
      (
        "number chunks",
        "# of chunks the file should be split into (conceptually), used to divide the file into disparate datasets",
        io::Serialization::GetAgent( &m_NumberChunks),
        "1"
      );

      member_data.AddInitializer
      (
        "chunks",
        "ranges of chunks to load, e.g. chunks=\"[ 0, 5) (7,10)\"",
        io::Serialization::GetAgent( &m_ChunkRanges),
        "[0]"
      );

      return member_data;
    } // GetParameters

    //! @brief get the size of the complete dataset without chunking
    //! @return the size of the complete dataset without chunking
    size_t RetrieveDatasetSubset::GetTotalDatasetSize() const
    {
      std::ifstream input;
      input.open( m_Filename.c_str(), std::ios::binary | std::ios::in);

      size_t byt_per_feature_result( 0);
      // skip over the feature and result labels
      {
        byt_per_feature_result += sizeof( float) * ReadFeatureLabelSetBinary( input).GetSize(); // feature
        byt_per_feature_result += sizeof( float) * ReadFeatureLabelSetBinary( input).GetSize(); // result
        byt_per_feature_result += ReadFeatureLabelSetBinary( input).GetSize(); // ids
      }

      std::ifstream::pos_type bytes_per_feature_result( byt_per_feature_result);

      // get the position before reading this line
      std::ifstream::pos_type first_feature_position( input.tellg());

      // determine the position at the end of the file
      input.seekg( 0, std::ios::end);

      std::ifstream::pos_type end_position( input.tellg());

      // close the file
      io::File::CloseClearFStream( input);

      // determine the # of chars in the file that are dedicated to the features/results
      std::ifstream::pos_type features_results_chars( end_position - first_feature_position);

      // Check that the size of the file for features is a multiple of line_width
      // Otherwise, the file is of the wrong format
      BCL_Assert
      (
        features_results_chars % bytes_per_feature_result == 0
        && first_feature_position > 0
        && bytes_per_feature_result > 0,
        "File " + m_Filename + " is of the wrong format; "
        + GetClassIdentifier() + " requires files to be in binary text format"
      );

      // determine the # of features
      return features_results_chars / bytes_per_feature_result;
    }

    //! @brief load a range of data from the dataset
    //! @param SUBSET the range of data to load
    //! @param FEATURES_STORAGE where to store features that are loaded, must be large enough to hold the subset without resizing
    //! @param RESULTS_STORAGE where to store the corresponding results, must be large enough to hold the subset without resizing
    //! @param START_FEATURE_NUMBER position to store the first feature in FEATURES_STORAGE
    //! @return # of features actually loaded
    //! Note: Implementations should overload this and SupportsEfficientSubsetLoading together
    size_t RetrieveDatasetSubset::GenerateDataSubset
    (
      const math::Range< size_t> &SUBSET,
      linal::MatrixInterface< float> &FEATURES_STORAGE,
      linal::MatrixInterface< float> &RESULTS_STORAGE,
      linal::MatrixInterface< char> &IDS_STORAGE,
      const size_t &START_FEATURE_NUMBER
    )
    {
      // read the data set from the file
      std::ifstream input;

      // open the file, determine the size, setup the feature and result columns
      const size_t total_dataset_size( GetTotalSizeAndParseLabels( input));

      // determine feature + result size
      const size_t feature_result_size( m_Features.GetInputFeatureSize());

      // determine id size
      const size_t ids_size( m_Ids.GetInputFeatureSize());

      // create a linal::Vector< float> object to hold the data loaded in
      linal::Vector< float> temp_vector( feature_result_size);
      linal::Vector< char> temp_id_vector( ids_size);

      // determine # of bytes per feature
      std::istream::pos_type bytes_per_feature( feature_result_size * sizeof( float) + ids_size);

      // get the actual ranges that will be taken from the dataset
      storage::Vector< math::Range< size_t> > ranges;
      ranges.AllocateMemory( m_ChunkRanges.GetRanges().GetSize());

      for
      (
        storage::Set< math::Range< size_t> >::const_iterator
          itr_range( m_ChunkRanges.GetRanges().Begin()), itr_range_end( m_ChunkRanges.GetRanges().End());
        itr_range != itr_range_end;
        ++itr_range
      )
      {
        ranges.PushBack
        (
          RetrieveDataSetBase::GetStartEndPositionOfRange( *itr_range, m_NumberChunks, total_dataset_size)
        );
      }

      // get the position before reading this line
      std::ifstream::pos_type first_feature_position( input.tellg());

      // keep track of how many rows have been loaded
      size_t loaded_feature_counter( 0), feature_position( START_FEATURE_NUMBER);

      // walk over the lines desired, loading only those columns and results that we want
      for
      (
        storage::Vector< math::Range< size_t> >::const_iterator itr_range( ranges.Begin()), itr_range_end( ranges.End());
        itr_range != itr_range_end && loaded_feature_counter < SUBSET.GetMax();
        ++itr_range
      )
      {
        if( itr_range->GetWidth() + loaded_feature_counter < SUBSET.GetMin())
        {
          loaded_feature_counter += itr_range->GetWidth();
          continue;
        }

        // determine 1st feature to load of this subset
        const size_t feature_offset
        (
          loaded_feature_counter >= SUBSET.GetMin()
          ? 0
          : SUBSET.GetMin() - loaded_feature_counter
        );

        loaded_feature_counter += feature_offset;

        // move to the appropriate line
        input.seekg
        (
          first_feature_position + std::istream::pos_type( itr_range->GetMin() + feature_offset) * bytes_per_feature
        );

        BCL_MessageVrb( "stream position: " + util::Format()( input.tellg()));

        // load in the desired lines
        for
        (
          size_t line_number( feature_offset), number_lines( itr_range->GetWidth());
          line_number < number_lines && loaded_feature_counter < SUBSET.GetMax();
          ++line_number, ++loaded_feature_counter
        )
        {
          // read the id vector
          io::BinarySerialize::ReadVector( temp_id_vector, input);
          // select the desired id columns
          m_Ids( temp_id_vector, IDS_STORAGE[ feature_position]);

          // read the features & results vector
          io::BinarySerialize::ReadVector( temp_vector, input);

          // put the chosen columns into the appropriate matrices
          m_Features( temp_vector, FEATURES_STORAGE[ feature_position]);
          m_Results( temp_vector, RESULTS_STORAGE[ feature_position]);
          ++feature_position;
        }
      }

      io::File::CloseClearFStream( input);

      return feature_position - START_FEATURE_NUMBER;
    }

    //! @brief get the nominal (e.g. best estimate without generating the entire dataset) size of the dataset
    //! @return the nominal (e.g. best estimate without generating the entire dataset) size of the dataset
    size_t RetrieveDatasetSubset::GetNominalSize() const
    {
      const size_t total_dataset_size( GetTotalDatasetSize());
      return RetrieveDataSetBase::GetSubsetSize( m_ChunkRanges, m_NumberChunks, total_dataset_size);
    }

    //! @brief read a feature label set written to a binary file
    FeatureLabelSet RetrieveDatasetSubset::ReadFeatureLabelSetBinary( std::ifstream &INPUT)
    {
      // get the position before reading this line
      std::ifstream::pos_type original_position( INPUT.tellg());

      FeatureLabelSet features;

      // read in the first string
      std::string class_name;
      INPUT >> class_name;
      if( class_name != features.GetClassIdentifier())
      {
        // bad string; probably an older format that had fewer feature label sets
        // reset the position
        INPUT.seekg( original_position);
        return features;
      }

      features.Read( INPUT);

      // read until finding the delimiter (null character).  This allows us to skip extraneous spaces
      if( INPUT.good() && !INPUT.eof())
      {
        char input( ' ');
        INPUT.get( input);
        while( !INPUT.eof() && input != '\0')
        {
          INPUT.get( input);
        }
      }

      // return the feature label
      return features;
    }

    //! @brief read a feature label set written to a binary file
    void RetrieveDatasetSubset::WriteFeatureLabelSetBinary( std::ostream &WRITE, const FeatureLabelSet &LABELS)
    {
      // write the string to output
      io::Serialize::Write( LABELS, WRITE);
      WRITE << '\0';
    }

    //! @brief read the feature and result labels and create DataSetSelectColumns from them
    //! @param INPUT the input stream to use; should not be open initially
    //!        after this function is called, it will be on the first line of the dataset
    //! @return the # of features in this file, which is determined here
    size_t RetrieveDatasetSubset::GetTotalSizeAndParseLabels( std::ifstream &INPUT)
    {
      io::File::CloseClearFStream( INPUT);
      INPUT.open( m_Filename.c_str(), std::ios::binary | std::ios::in);

      // load the feature and result labels
      FeatureLabelSet feature_labels, result_labels, id_labels;
      feature_labels = ReadFeatureLabelSetBinary( INPUT);
      result_labels = ReadFeatureLabelSetBinary( INPUT);
      id_labels = ReadFeatureLabelSetBinary( INPUT);

      // determine the total # of features/results on each line
      const size_t total_features_results( feature_labels.GetSize() + result_labels.GetSize());

      // determine # of features
      size_t number_features( 0);

      {
        // get the position of the first feature-result
        std::istream::pos_type first_feature_position( INPUT.tellg());

        // compute bytes per feautre
        std::istream::pos_type bytes_per_feature_result( sizeof( float) * total_features_results + id_labels.GetSize());

        // determine the position at the end of the file
        INPUT.seekg( 0, std::ios::end);
        std::ifstream::pos_type end_position( INPUT.tellg());

        // return the stream to the first feature
        INPUT.seekg( first_feature_position);

        // determine the # of bytes in the file that are dedicated to the features/results
        std::ifstream::pos_type features_results_bytes( end_position - first_feature_position);

        BCL_MessageDbg( "feature_labels: " + util::Format()( feature_labels));

        BCL_MessageDbg( "result_labels: " + util::Format()( result_labels));

        BCL_MessageDbg
        (
          "features_results_bytes: " + util::Format()( features_results_bytes) +
          "\nbytes_per_feature_result: " + util::Format()( bytes_per_feature_result) +
          "\nfirst_feature_position: " + util::Format()( first_feature_position) +
          "\nbytes_per_feature_result: " + util::Format()( bytes_per_feature_result)
        );

        // Check that the size of the file for features is a multiple of line_width
        // Otherwise, the file is of the wrong format
        BCL_Assert
        (
          features_results_bytes % bytes_per_feature_result == 0
          && first_feature_position > 0
          && bytes_per_feature_result > 0,
          "File " + m_Filename + " is of the wrong format; "
          + GetClassIdentifier() + " requires files to be in binary text format"
        );

        // determine the # of features
        number_features = features_results_bytes / bytes_per_feature_result;
      }

      // check if the features selector has not already been set up
      if( !m_Features.GetInputFeatureSize())
      {
        // determine which feature labels will be kept
        storage::Vector< size_t> features_kept;
        features_kept.AllocateMemory( feature_labels.GetSize());

        // check; if no features were given, then set all features
        if( GetFeatureCode().IsScalar())
        {
          features_kept.AllocateMemory( feature_labels.GetSize());
          for( size_t i( 0), feature_size( feature_labels.GetSize()); i < feature_size; ++i)
          {
            features_kept.PushBack( i);
          }
        }
        else
        {
          // walk over the arguments in the feature label
          for
          (
            util::ObjectDataLabel::const_iterator
              itr_label( GetFeatureCode().Begin()), itr_label_end( GetFeatureCode().End());
            itr_label != itr_label_end;
            ++itr_label
          )
          {
            // append the indices of the given feature
            features_kept.Append( feature_labels.GetPropertyIndices( *itr_label));
          }
        }
        // create an appropriate DataSetSelectColumns object to select the desired columns of the dataset
        m_Features = DataSetSelectColumns( total_features_results, features_kept);
        BCL_MessageDbg( "Feature selector: " + util::Format()( m_Features));
      }

      if( !m_Results.GetInputFeatureSize())
      {
        // determine which feature labels will be kept
        storage::Vector< size_t> results_kept;
        results_kept.AllocateMemory( result_labels.GetSize());

        // check; if no result code given, then set all results
        if( GetResultCode().IsScalar())
        {
          results_kept.AllocateMemory( result_labels.GetSize());
          for( size_t i( 0), result_size( result_labels.GetSize()); i < result_size; ++i)
          {
            results_kept.PushBack( i);
          }
        }
        else
        {
          // now add the results labels, remembering to always add the size of the features first
          for
          (
            util::ObjectDataLabel::const_iterator
              itr_label( GetResultCode().Begin()), itr_label_end( GetResultCode().End());
            itr_label != itr_label_end;
            ++itr_label
          )
          {
            // append the indices of the given result
            results_kept.Append( result_labels.GetPropertyIndices( *itr_label));
          }
        }
        // add the offset to the results kept vector, since results come after features
        for
        (
          storage::Vector< size_t>::iterator itr( results_kept.Begin()), itr_end( results_kept.End());
          itr != itr_end;
          ++itr
        )
        {
          *itr += feature_labels.GetSize();
        }

        m_Results = DataSetSelectColumns( total_features_results, results_kept);
        BCL_MessageDbg( "Result selector: " + util::Format()( m_Results));
      }

      if( !m_Ids.GetInputFeatureSize())
      {
        // determine which feature labels will be kept
        storage::Vector< size_t> ids_kept;
        ids_kept.AllocateMemory( id_labels.GetSize());

        // check; if no result code given, then set all results
        if( GetIdCode().IsScalar())
        {
          ids_kept.AllocateMemory( id_labels.GetSize());
          for( size_t i( 0), id_size( id_labels.GetSize()); i < id_size; ++i)
          {
            ids_kept.PushBack( i);
          }
        }
        else
        {
          // now add the results labels, remembering to always add the size of the features first
          for
          (
            util::ObjectDataLabel::const_iterator
              itr_label( GetIdCode().Begin()), itr_label_end( GetIdCode().End());
            itr_label != itr_label_end;
            ++itr_label
          )
          {
            // append the indices of the given result
            ids_kept.Append( id_labels.GetPropertyIndices( *itr_label));
          }
        }

        m_Ids = DataSetSelectColumns( id_labels.GetSize(), ids_kept);
        BCL_MessageDbg( "IDs selector: " + util::Format()( m_Ids));
      }
      BCL_MessageStd
      (
        "# features: " + util::Format()( number_features)
        + " feature size: " + util::Format()( m_Features.GetOutputFeatureSize())
        + " result size: " + util::Format()( m_Results.GetOutputFeatureSize())
        + " feature result size of data superset: " + util::Format()( m_Results.GetInputFeatureSize())
      );
      return number_features;
    }

  } // namespace model
} // namespace bcl
