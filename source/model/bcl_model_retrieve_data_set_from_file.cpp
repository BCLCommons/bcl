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
#include "model/bcl_model_retrieve_data_set_from_file.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_file.h"
#include "io/bcl_io_serialization.h"
#include "linal/bcl_linal_vector_operations.h"
#include "math/bcl_math_statistics.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {

  //////////
  // data //
  //////////

    // add the interface to the set of known implementations
    const util::SiPtr< const util::ObjectInterface> RetrieveDataSetFromFile::s_Instance
    (
      util::Enumerated< RetrieveDataSetBase>::AddInstance( new RetrieveDataSetFromFile())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    //! @param FILENAME the file to retrieve the data set from
    //! @param NUMBER_CHUNKS # of chunks the file should be split into (conceptually), used to divide the file into disparate datasets
    //! @param CHUNKS chunks to load
    RetrieveDataSetFromFile::RetrieveDataSetFromFile
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
    RetrieveDataSetFromFile *RetrieveDataSetFromFile::Clone() const
    {
      return new RetrieveDataSetFromFile( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &RetrieveDataSetFromFile::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &RetrieveDataSetFromFile::GetAlias() const
    {
      static const std::string s_Name( "File");
      return s_Name;
    }

    //! @brief Get the code / label set for the feature (1st part) of the data set with sizes of each property
    //! @return the code / label for the feature (1st part) of the data set with sizes of each property
    //! the feature code set
    FeatureLabelSet RetrieveDataSetFromFile::GetFeatureLabelsWithSizes() const
    {
      storage::Vector< util::ObjectDataLabel> labels;

      const size_t feature_size( GetFeatureResultSizes().First());

      // create 1 feature label with size 1 for each feature
      storage::Vector< size_t> sizes( feature_size, 1);

      if( GetFeatureCode().GetNumberArguments() == size_t( 0))
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
        sizes.Resize( labels.GetSize());
      }

      return FeatureLabelSet( "Combine", labels, sizes);
    }

    //! @brief Get the code / label for the result (2nd part) of the data set with sizes of each property
    //! @return the code / label for the result (2nd part) of the data set with sizes of each property
    //! the feature code set
    FeatureLabelSet RetrieveDataSetFromFile::GetResultCodeWithSizes() const
    {
      storage::Vector< util::ObjectDataLabel> labels;
      const storage::VectorND< 2, size_t> feature_result_sizes( GetFeatureResultSizes());
      const size_t original_feature_columns
      (
        m_Features.GetInputFeatureSize()
        ? m_Features.GetInputFeatureSize()
        : feature_result_sizes.First()
      );
      const size_t result_size( feature_result_sizes.Second());

      // create 1 label with size 1 for each result
      storage::Vector< size_t> sizes( result_size, 1);

      if( GetResultCode().GetNumberArguments() == size_t( 0))
      {
        // the labels will just be the index of each feature
        labels.AllocateMemory( result_size);
        for( size_t result_number( 0); result_number < result_size; ++result_number)
        {
          labels.PushBack( util::ObjectDataLabel( util::Format()( original_feature_columns + result_number)));
        }
      }
      else
      {
        // use the given labels
        labels = GetResultCode().GetArguments();
        sizes.Resize( labels.GetSize());
      }

      return FeatureLabelSet( "Combine", labels, sizes);
    }

    //! @brief Set the code / label for the feature (1st part) of the data set
    //! @param CODE the new code
    void RetrieveDataSetFromFile::SelectFeatures( const util::ObjectDataLabel &CODE)
    {
      const size_t original_feature_columns
      (
        m_Features.GetInputFeatureSize()
        ? m_Features.GetInputFeatureSize()
        : GetFeatureResultSizes().First()
      );
      if( CODE.GetNumberArguments() == original_feature_columns)
      {
        return;
      }
      storage::Vector< size_t>
        indices( util::SplitStringToNumerical< size_t>( CODE.ArgumentsToString( false, ','), ","));
      BCL_Assert
      (
        math::Statistics::MaximumValue( indices.Begin(), indices.End()) < original_feature_columns,
        "SelectFeatures given index beyond end of array"
      );
      m_Features = DataSetSelectColumns( original_feature_columns, indices);
      RetrieveDataSetBase::SelectFeatures( CODE);
    }

    //! @brief Set the code / label for the result (2nd part) of the data set
    //! @param CODE the new code
    void RetrieveDataSetFromFile::SelectResults( const util::ObjectDataLabel &CODE)
    {
      const size_t original_feature_columns
      (
        m_Features.GetInputFeatureSize()
        ? m_Features.GetInputFeatureSize()
        : GetFeatureResultSizes().First()
      );
      const size_t original_result_columns
      (
        m_Results.GetInputFeatureSize()
        ? m_Results.GetInputFeatureSize()
        : GetFeatureResultSizes().Second()
      );
      if( CODE.GetNumberArguments() == original_result_columns)
      {
        return;
      }
      storage::Vector< size_t>
        indices( util::SplitStringToNumerical< size_t>( CODE.ArgumentsToString( false, ','), ","));
      BCL_Assert
      (
        math::Statistics::MaximumValue( indices.Begin(), indices.End())
          < original_feature_columns + original_result_columns,
        "SelectResults given index beyond end of array"
      );
      BCL_Assert
      (
        math::Statistics::MinimumValue( indices.Begin(), indices.End()) >= original_feature_columns,
        "SelectResults given feature index"
      );
      linal::VectorReference< size_t> indices_ref( indices.GetSize(), &*indices.Begin());
      indices_ref -= original_feature_columns;
      m_Features = DataSetSelectColumns( original_result_columns, indices);
      RetrieveDataSetBase::SelectResults( CODE);
    }

    //! @brief Get the code / label for the ids of the data set with sizes of each property
    //! @return the code / label for the ids of the data set with sizes of each property
    //! the feature code set
    FeatureLabelSet RetrieveDataSetFromFile::GetIdCodeWithSizes() const
    {
      // files currently have no way of storing ids
      return FeatureLabelSet();
    }

    //! @brief get the number of partitions requested by the user, along with the partition ids
    //! @return the number of partitions requested by the user, along with the partition ids
    storage::Pair< size_t, math::RangeSet< size_t> > RetrieveDataSetFromFile::GetNumberPartitionsAndIds() const
    {
      return storage::Pair< size_t, math::RangeSet< size_t> >( m_NumberChunks, m_ChunkRanges);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief generate dataset
    //! @return generated dataset
    util::ShPtr< descriptor::Dataset>
    RetrieveDataSetFromFile::GenerateDataSet()
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

      // determine the # of feature/results
      const size_t total_number_feature_results( GetTotalDatasetSize());

      BCL_MessageDbg( "file: " + m_Filename + " total # features: " + util::Format()( total_number_feature_results));

      // determine the # of features/results in the set of data that will be loaded
      const size_t number_feature_results
      (
        RetrieveDataSetBase::GetSubsetSize( m_ChunkRanges, m_NumberChunks, total_number_feature_results)
      );

      // create a vector that will store bools that indicate whether or not to load each row
      storage::Vector< size_t> should_load_row( total_number_feature_results, size_t( 0));

      // set the vector up with the given ranges
      for
      (
        storage::Set< math::Range< size_t> >::const_iterator
          itr( m_ChunkRanges.GetRanges().Begin()), itr_end( m_ChunkRanges.GetRanges().End());
        itr != itr_end;
        ++itr
      )
      {
        math::Range< size_t> actual_range
        (
          RetrieveDataSetBase::GetStartEndPositionOfRange( *itr, m_NumberChunks, total_number_feature_results)
        );
        for
        (
          size_t row_id( actual_range.GetMin()), max_row_id( actual_range.GetMax());
          row_id < max_row_id;
          ++row_id
        )
        {
          should_load_row( row_id) = 1;
        }
      }

      // open up the
      // read the data set from the file
      io::IFStream input;
      io::File::MustOpenIFStream( input, m_Filename);

      // look for this string; the number that follows is the total size
      const std::string search_string
      (
        GetStaticClassName< storage::VectorND< 2, linal::Vector< double> > >()
      );

      storage::VectorND< 2, size_t> feature_result_sizes( GetFeatureResultSizes());

      // search for search string, which will be the first line of the vector
      std::string temp_string;
      while( input.good() && !input.eof())
      {
        std::streampos old_position( input.tellg());
        std::getline( input, temp_string);
        if( util::TrimString( temp_string) == search_string)
        {
          // move back to the position right before the search string
          input.seekg( old_position);
          break;
        }
      }

      // allocate enough memory for the features that will be loaded
      util::ShPtr< descriptor::Dataset> dataset
      (
        new descriptor::Dataset
        (
          number_feature_results,
          GetFeatureLabelsWithSizes(),
          GetResultCodeWithSizes(),
          GetIdCodeWithSizes()
        )
      );

      linal::MatrixReference< float> features( dataset->GetFeaturesReference());
      linal::MatrixReference< float> results( dataset->GetResultsReference());

      // read in the desired rows
      storage::VectorND< 2, linal::Vector< double> > temp_row;
      linal::Vector< float>
        temp_feature( m_Features.GetInputFeatureSize()), temp_result( m_Results.GetInputFeatureSize());
      // the results need to go through DataSetSelectColumns
      for( size_t row_number( 0), feature_number( 0); row_number < total_number_feature_results; ++row_number)
      {
        input >> temp_row;
        if( should_load_row( row_number))
        {
          // filter the features, if only selecting some columns
          if( m_Features.GetInputFeatureSize())
          {
            // filter the result
            std::copy( temp_row.First().Begin(), temp_row.First().End(), temp_feature.Begin());
            m_Features( temp_feature, features[ feature_number]);
          }
          else
          {
            std::copy( temp_row.First().Begin(), temp_row.First().End(), features[ feature_number]);
          }

          // filter the result, if only selecting some columns
          if( m_Results.GetInputFeatureSize())
          {
            // filter the result
            std::copy( temp_row.Second().Begin(), temp_row.Second().End(), temp_result.Begin());
            m_Results( temp_result, results[ feature_number]);
          }
          else
          {
            // copy all columns of the result
            std::copy( temp_row.Second().Begin(), temp_row.Second().End(), results[ feature_number]);
          }
          ++feature_number;
        }
      }

      io::File::CloseClearFStream( input);

      // return the generated data set
      return dataset;
    } // GenerateDataSet

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer RetrieveDataSetFromFile::GetSerializer() const
    {
      io::Serializer member_data;
      member_data.SetClassDescription
      (
        "retrieves a data set from a file containing a single " +
        GetStaticClassName< storage::Vector< storage::VectorND< 2, linal::Vector< double> > > >()
      );

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
    size_t RetrieveDataSetFromFile::GetTotalDatasetSize() const
    {
      // read the data set from the file
      io::IFStream input;
      io::File::MustOpenIFStream( input, m_Filename);

      // look for this string; the number that follows is the total size
      const std::string search_string
      (
        GetStaticClassName< storage::Vector< storage::VectorND< 2, linal::Vector< double> > > >()
      );

      // search for search string, which will be the first line of the vector
      std::string temp_string;
      while( input.good() && !input.eof())
      {
        std::getline( input, temp_string);
        if( util::TrimString( temp_string) == search_string)
        {
          break;
        }
      }

      // the number that follows is the total size
      size_t size( 0);
      if( input.good() && !input.eof())
      {
        input >> size;
      }
      else
      {
        BCL_MessageCrt
        (
          "Dataset will be empty, could not find the expected dataset class identifier: " + search_string
        );
      }

      io::File::CloseClearFStream( input);

      return size;
    }

    //! @brief get the nominal (e.g. best estimate without generating the entire dataset) size of the dataset
    //! @return the nominal (e.g. best estimate without generating the entire dataset) size of the dataset
    size_t RetrieveDataSetFromFile::GetNominalSize() const
    {
      const size_t total_dataset_size( GetTotalDatasetSize());
      return RetrieveDataSetBase::GetSubsetSize( m_ChunkRanges, m_NumberChunks, total_dataset_size);
    }

    //! @brief get the sizes of the features and results in the dataset without loading the whole dataset
    //! @return the sizes of the features and results in the dataset without loading the whole dataset
    storage::VectorND< 2, size_t> RetrieveDataSetFromFile::GetFeatureResultSizes() const
    {
      // first, check if only selected columns have been enabled; if so, take those sizes rather than the complete
      // size of the dataset
      const size_t columns_size_feature( m_Features.GetOutputFeatureSize());
      const size_t columns_size_results( m_Results.GetOutputFeatureSize());

      if( m_Features.GetInputFeatureSize() && m_Results.GetInputFeatureSize())
      {
        // both the dataset select columns were enabled, return the selected sizes
        return storage::VectorND< 2, size_t>( columns_size_feature, columns_size_results);
      }

      // read the data set from the file
      io::IFStream input;
      io::File::MustOpenIFStream( input, m_Filename);

      // look for this string; the number that follows is the nominal size
      const std::string search_string
      (
        GetStaticClassName< storage::VectorND< 2, linal::Vector< double> > >()
      );

      // search for search string, which will be the first line of the vector
      std::string temp_string;
      while( input.good() && !input.eof())
      {
        std::getline( input, temp_string);
        if( util::TrimString( temp_string) == search_string)
        {
          break;
        }
      }

      storage::VectorND< 2, size_t> feature_result_sizes
      (
        util::GetUndefined< size_t>(),
        util::GetUndefined< size_t>()
      );

      // the vectors that follow contain the first features and results
      if( input.good() && !input.eof())
      {
        linal::Vector< double> temp_vector;
        input >> temp_vector;
        feature_result_sizes.First() = temp_vector.GetSize();
        input >> temp_vector;
        feature_result_sizes.Second() = temp_vector.GetSize();
      }

      // override the given sizes by the sizes of the selects, if available
      if( m_Features.GetInputFeatureSize())
      {
        feature_result_sizes.First() = columns_size_feature;
      }

      if( m_Results.GetInputFeatureSize())
      {
        feature_result_sizes.Second() = columns_size_results;
      }

      io::File::CloseClearFStream( input);

      return feature_result_sizes;
    }

  } // namespace model
} // namespace bcl
