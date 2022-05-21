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
#include "model/bcl_model_retrieve_data_set_from_delimited_file.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_file.h"
#include "io/bcl_io_serialization.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {

  //////////
  // data //
  //////////

    // add the interface to the set of known implementations
    const util::SiPtr< const util::ObjectInterface> RetrieveDataSetFromDelimitedFile::s_Instance
    (
      util::Enumerated< RetrieveDataSetBase>::AddInstance( new RetrieveDataSetFromDelimitedFile())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    //! @param FILENAME the file to retrieve the data set from
    //! @param NUMBER_CHUNKS # of chunks the file should be split into (conceptually), used to divide the file into disparate datasets
    //! @param CHUNKS chunks to load
    RetrieveDataSetFromDelimitedFile::RetrieveDataSetFromDelimitedFile
    (
      const std::string &FILENAME,
      const size_t &NUMBER_RESULT_COLS,
      const size_t &NUMBER_ID_CHARS
    ) :
      m_Filename( FILENAME),
      m_Features(),
      m_Results(),
      m_Ids(),
      m_ResultsByNumberLastColumns( NUMBER_RESULT_COLS),
      m_NumberIdCharacters( NUMBER_ID_CHARS)
    {
    }

    //! @brief Clone is the virtual Clone constructor
    //! @return a pointer to new LocatorAtomFurthest which is a copy of this
    RetrieveDataSetFromDelimitedFile *RetrieveDataSetFromDelimitedFile::Clone() const
    {
      return new RetrieveDataSetFromDelimitedFile( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &RetrieveDataSetFromDelimitedFile::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &RetrieveDataSetFromDelimitedFile::GetAlias() const
    {
      static const std::string s_Name( "Csv");
      return s_Name;
    }

    //! @brief Get the code / label set for the feature (1st part) of the data set with sizes of each property
    //! @return the code / label for the feature (1st part) of the data set with sizes of each property
    //! the feature code set
    FeatureLabelSet RetrieveDataSetFromDelimitedFile::GetFeatureLabelsWithSizes() const
    {
      storage::Vector< util::ObjectDataLabel> labels;

      const size_t feature_size( GetFeatureResultSizes().First());

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
    FeatureLabelSet RetrieveDataSetFromDelimitedFile::GetResultCodeWithSizes() const
    {
      storage::Vector< util::ObjectDataLabel> labels;

      const size_t result_size( GetFeatureResultSizes().Second());

      // create 1 label with size 1 for each result
      storage::Vector< size_t> sizes( result_size, 1);

      if( GetResultCode().GetNumberArguments() != result_size)
      {
        // the labels will just be the index of each feature
        labels.AllocateMemory( result_size);
        for( size_t result_number( 0); result_number < result_size; ++result_number)
        {
          labels.PushBack( util::ObjectDataLabel( util::Format()( result_number)));
        }
      }
      else
      {
        // use the given labels
        labels = GetResultCode().GetArguments();
      }

      return FeatureLabelSet( "Combine", labels, sizes);
    }

    //! @brief Get the code / label for the ids of the data set with sizes of each property
    //! @return the code / label for the ids of the data set with sizes of each property
    //! the feature code set
    FeatureLabelSet RetrieveDataSetFromDelimitedFile::GetIdCodeWithSizes() const
    {
      if( !m_NumberIdCharacters)
      {
        return FeatureLabelSet();
      }
      storage::Vector< util::ObjectDataLabel> labels;

      // create 1 label with size 1 for each result
      storage::Vector< size_t> sizes( m_NumberIdCharacters, size_t( 1));

      if( GetIdCode().GetNumberArguments() != m_NumberIdCharacters)
      {
        // the labels will just be the index of each feature
        labels.AllocateMemory( m_NumberIdCharacters);
        for( size_t id_number( 0); id_number < m_NumberIdCharacters; ++id_number)
        {
          labels.PushBack( util::ObjectDataLabel( util::Format()( id_number)));
        }
      }
      else
      {
        // use the given labels
        labels = GetIdCode().GetArguments();
      }

      return FeatureLabelSet( "Combine", labels, sizes);
    }

    //! @brief get the number of partitions requested by the user, along with the partition ids
    //! @return the number of partitions requested by the user, along with the partition ids
    storage::Pair< size_t, math::RangeSet< size_t> > RetrieveDataSetFromDelimitedFile::GetNumberPartitionsAndIds() const
    {
      return storage::Pair< size_t, math::RangeSet< size_t> >( 1, math::RangeSet< size_t>( math::Range< size_t>( 0, 0)));
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief generate dataset
    //! @return generated dataset
    util::ShPtr< descriptor::Dataset>
    RetrieveDataSetFromDelimitedFile::GenerateDataSet()
    {
      // determine the # of feature/results
      const size_t total_number_feature_results( GetTotalDatasetSize());

      static const std::string s_delimiters( " ,\t|");

      BCL_MessageDbg( "file: " + m_Filename + " total # features: " + util::Format()( total_number_feature_results));

      // open up the
      // read the data set from the file
      io::IFStream input;
      io::File::MustOpenIFStream( input, m_Filename);

      // allocate enough memory for the features that will be loaded
      util::ShPtr< descriptor::Dataset> dataset
      (
        new descriptor::Dataset
        (
          total_number_feature_results,
          GetFeatureLabelsWithSizes(),
          GetResultCodeWithSizes(),
          GetIdCodeWithSizes()
        )
      );

      linal::MatrixReference< float> features( dataset->GetFeaturesReference());
      linal::MatrixReference< float> results( dataset->GetResultsReference());
      linal::MatrixReference< char> ids( dataset->GetIdsReference());

      // read in the desired rows
      storage::VectorND< 2, linal::Vector< double> > temp_row;
      linal::Vector< float>
        temp_feature( features.GetNumberCols()), temp_result( results.GetNumberCols());

      // the results need to go through DataSetSelectColumns

      std::string temp_line;

      storage::Vector< float> temp_line_vector;

      const size_t results_position( temp_feature.GetSize());

      for( size_t row_number( 0); row_number < total_number_feature_results; ++row_number)
      {
        std::getline( input, temp_line);

        if( m_NumberIdCharacters)
        {
          // prune off the id
          std::copy( temp_line.begin(), temp_line.begin() + m_NumberIdCharacters, ids[ row_number]);
          BCL_Assert
          (
            s_delimiters.find( temp_line[ m_NumberIdCharacters]) != std::string::npos,
            "IDs should be fixed width, followed by a valid delimiter! Failed on line: " + util::Format()( row_number)
            + " with line: " + temp_line
          );
          temp_line = temp_line.substr( m_NumberIdCharacters + 1);
        }

        temp_line_vector = util::SplitStringToNumerical< float>( temp_line, s_delimiters);

        // determine
        storage::Vector< float>::iterator features_end_itr( temp_line_vector.Begin() + results_position);
        BCL_Assert
        (
          temp_line_vector.GetSize()
          ==
          m_Results.GetInputFeatureSize() + m_Features.GetInputFeatureSize()
          + features.GetNumberCols() + results.GetNumberCols(),
          "Line #" + util::Format()( row_number) + " (0-indexed) with incorrect number of values; had: "
          + util::Format()( temp_line_vector.GetSize()) + " but expected: " +
          util::Format()
          (
              m_Results.GetInputFeatureSize() + m_Features.GetInputFeatureSize()
              + features.GetNumberCols() + results.GetNumberCols()
           )
        );

        // filter the features, if only selecting some columns
        if( m_Features.GetInputFeatureSize())
        {
          // filter the result
          std::copy( temp_line_vector.Begin(), features_end_itr, temp_feature.Begin());
          m_Features( temp_feature, features[ row_number]);
        }
        else
        {
          std::copy( temp_line_vector.Begin(), features_end_itr, features[ row_number]);
        }

        // filter the result, if only selecting some columns
        if( m_Results.GetInputFeatureSize())
        {
          // filter the result
          std::copy( features_end_itr, temp_line_vector.End(), temp_result.Begin());
          m_Results( temp_result, results[ row_number]);
        }
        else
        {
          // copy all columns of the result
          std::copy( features_end_itr, temp_line_vector.End(), results[ row_number]);
        }
      }

      io::File::CloseClearFStream( input);

      // return the generated data set
      return dataset;
    } // GenerateDataSet

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer RetrieveDataSetFromDelimitedFile::GetSerializer() const
    {
      io::Serializer member_data;
      member_data.SetClassDescription
      (
        "retrieves a data set from a delimited value file. Possible delimiters are comma, tab, space, and pipe (|)"
      );

      member_data.AddInitializer
      (
        "filename",
        "name of the file from which to load the dataset",
        io::Serialization::GetAgentInputFilename( &m_Filename)
      );

      member_data.AddInitializer
      (
        "number result cols",
        "number of columns (to the right of the data set) that specify the results",
        io::Serialization::GetAgent( &m_ResultsByNumberLastColumns),
        "1"
      );

      member_data.AddInitializer
      (
        "number id chars",
        "number of characters (to the left of the data set) that specify the ids. "
        "There should be a delimiter (which is ignored) between the ids and features, if any ids are present",
        io::Serialization::GetAgent( &m_NumberIdCharacters),
        "0"
      );

      return member_data;
    } // GetParameters

    //! @brief get the size of the complete dataset without chunking
    //! @return the size of the complete dataset without chunking
    size_t RetrieveDataSetFromDelimitedFile::GetTotalDatasetSize() const
    {
      // read the data set from the file
      io::IFStream input;
      io::File::MustOpenIFStream( input, m_Filename);

      // search for search string, which will be the first line of the vector
      size_t counter( 0);
      std::string temp_string;
      while( input.good())
      {
        std::getline( input, temp_string);

        // if eof check is in while condition it will count eof as an additional line
        if( !input.eof())
        {
          ++counter;
        }
        else
        {
          break;
        }
      }

      // if there where no lines in your delimited data file
      if( counter == 0)
      {
        BCL_MessageStd( "No rows available in delimited data file!");
      }

      io::File::CloseClearFStream( input);

      return counter;
    }

    //! @brief get the nominal (e.g. best estimate without generating the entire dataset) size of the dataset
    //! @return the nominal (e.g. best estimate without generating the entire dataset) size of the dataset
    size_t RetrieveDataSetFromDelimitedFile::GetNominalSize() const
    {
      return GetTotalDatasetSize();
    }

    //! @brief get the sizes of the features and results in the dataset without loading the whole dataset
    //! @return the sizes of the features and results in the dataset without loading the whole dataset
    storage::VectorND< 2, size_t> RetrieveDataSetFromDelimitedFile::GetFeatureResultSizes() const
    {
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

      storage::VectorND< 2, size_t> feature_result_sizes
      (
        util::GetUndefined< size_t>(),
        util::GetUndefined< size_t>()
      );

      // the vectors that follow contain the first features and results
      if( input.good() && !input.eof())
      {
        // extract line and tokenize
        std::string temp_vector;
        std::getline( input, temp_vector);
        if( !input.good() && util::TrimString( temp_vector).empty())
        {
          return storage::VectorND< 2, size_t>( columns_size_feature, columns_size_results);
        }
        if( m_NumberIdCharacters)
        {
          // prune off the id part, if any was present
          BCL_Assert( m_NumberIdCharacters < temp_vector.size() + 1, "More ids than characters on line: " + temp_vector);
          temp_vector = temp_vector.substr( m_NumberIdCharacters + 1);
        }
        storage::Vector< std::string> split_values( util::SplitString( temp_vector, " ,\t|"));

        BCL_Assert( m_ResultsByNumberLastColumns < split_values.GetSize(), "More results columns specified than columns available! Perhaps delimiter is wrong?");

        feature_result_sizes.First() = split_values.GetSize() - m_ResultsByNumberLastColumns;
        feature_result_sizes.Second() = m_ResultsByNumberLastColumns;
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

      // close file stream
      io::File::CloseClearFStream( input);

      return feature_result_sizes;
    }

  } // namespace model
} // namespace bcl
