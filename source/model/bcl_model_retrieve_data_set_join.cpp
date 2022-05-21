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
#include "model/bcl_model_retrieve_data_set_join.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_file.h"
#include "io/bcl_io_ifstream.h"
#include "io/bcl_io_serialization.h"
#include "storage/bcl_storage_triplet.h"
// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {

  //////////
  // data //
  //////////

    // add the interface to the set of known implementations
    const util::SiPtr< const util::ObjectInterface> RetrieveDataSetJoin::s_Instance
    (
      util::Enumerated< RetrieveDataSetBase>::AddInstance( new RetrieveDataSetJoin())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief Clone is the virtual Clone constructor
    //! @return a pointer to new LocatorAtomFurthest which is a copy of this
    RetrieveDataSetJoin *RetrieveDataSetJoin::Clone() const
    {
      return new RetrieveDataSetJoin( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &RetrieveDataSetJoin::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &RetrieveDataSetJoin::GetAlias() const
    {
      static const std::string s_Name( "Join");
      return s_Name;
    }

    void RetrieveDataSetJoin::SelectFeatures( const util::ObjectDataLabel &CODE)
    {
      RetrieveDataSetBase::SelectFeatures( CODE);
    }

    //! @brief Set the code / label for the result (2nd part) of the data set
    //! @return the code / label for the result (2nd part) of the data set
    void RetrieveDataSetJoin::SelectResults( const util::ObjectDataLabel &CODE)
    {
      RetrieveDataSetBase::SelectResults( CODE);
    }

    //! @brief Set which id columns to retrieve
    //! @param CODE the id column names to retrieve
    void RetrieveDataSetJoin::SelectIds( const util::ObjectDataLabel &CODE)
    {
      RetrieveDataSetBase::SelectIds( CODE);
//      m_LHS->SelectIds( CODE);
    }

    //! @brief Get the code / label set for the feature (1st part) of the data set with sizes of each property
    //! @return the code / label for the feature (1st part) of the data set with sizes of each property
    //! the feature code set
    FeatureLabelSet RetrieveDataSetJoin::GetFeatureLabelsWithSizes() const
    {
//      FeatureLabelSet label_set( m_LHS->GetFeatureLabelsWithSizes());
//      label_set.
    }

    //! @brief Get the code / label for the result (2nd part) of the data set with sizes of each property
    //! @return the code / label for the result (2nd part) of the data set with sizes of each property
    //! the feature code set
    FeatureLabelSet RetrieveDataSetJoin::GetResultCodeWithSizes() const
    {

//      return GetFeatureLabelsWithSizes();
    }

    //! @brief Get the code / label for the ids of the data set with sizes of each property
    //! @return the code / label for the ids of the data set with sizes of each property
    //! the feature code set
    FeatureLabelSet RetrieveDataSetJoin::GetIdCodeWithSizes() const
    {
      return m_LHS->GetIdCodeWithSizes();
    }

    //! @brief get the number of partitions requested by the user, along with the partition ids
    //! @return the number of partitions requested by the user, along with the partition ids
    storage::Pair< size_t, math::RangeSet< size_t> > RetrieveDataSetJoin::GetNumberPartitionsAndIds() const
    {
      return m_LHS->GetNumberPartitionsAndIds();
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief generate dataset
    //! @return generated dataset
    util::ShPtr< descriptor::Dataset>
    RetrieveDataSetJoin::GenerateDataSet()
    {
      UpdateIDsIndex();
      // determine the approximate size of the complete dataset after join
      const size_t full_dataset_size( GetNominalSize());

      // we will need to map id label pairs to results
      storage::Vector< storage::Triplet< std::string, std::string, linal::Vector< float> > > id_pairs_to_result;
      std::string temp_line;
      storage::Vector< std::string> temp_line_vector;
      storage::Pair< std::string, std::string> temp_pair;

      // read id labels corresponding to datasets to be joined
      io::IFStream input;
      io::File::MustOpenIFStream( input, m_IdPairsFileName);
      for( size_t id_row_number( 0); id_row_number < full_dataset_size; ++id_row_number)
      {
        // parse each line for datasets to be joined
        std::getline( input, temp_line);
        id_pairs_to_result.PushBack();
        temp_line_vector = util::SplitString( temp_line, ",\r\n");
        id_pairs_to_result.LastElement().First() = temp_line_vector( 0);
        id_pairs_to_result.LastElement().Second() = temp_line_vector( 1);

        linal::Vector< float> &temp_result_vector( id_pairs_to_result.LastElement().Third());
        temp_result_vector = linal::Vector< float>( temp_line_vector.GetSize() - 2);
        BCL_Assert
        (
          temp_result_vector.GetSize() == id_pairs_to_result.FirstElement().Third().GetSize(),
          "Different # of results for various id pairs on line #" + util::Format()( id_row_number + 1)
        );
        for( size_t i( 0), sz( temp_result_vector.GetSize()); i < sz; ++i)
        {
          temp_result_vector( i) = util::ConvertStringToNumericalValue< float>( temp_line_vector( i));
        }
      }
      const size_t result_size( id_pairs_to_result.LastElement().Third().GetSize());

      // close the input file stream
      io::File::CloseClearFStream( input);

      util::ShPtr< descriptor::Dataset> lhs_dataset_ptr( m_LHS->GenerateDataSet());
      util::ShPtr< descriptor::Dataset> rhs_dataset_ptr( m_RHS->GenerateDataSet());

      storage::Map< std::string, storage::Pair< std::string, size_t>> lhs_map_property_to_index, rhs_map_property_to_index;
      BCL_Assert( m_LHSIdIndices.LastElement() < lhs_dataset_ptr->GetIdsReference().GetNumberCols(), "LHS out of bounds id column!");
      BCL_Assert( m_RHSIdIndices.LastElement() < rhs_dataset_ptr->GetIdsReference().GetNumberCols(), "RHS out of bounds id column!");
      const size_t lhs_n_id_cols( m_LHSIdIndices.GetSize()), rhs_n_id_cols( m_RHSIdIndices.GetSize());
      std::string lhs_tmp( m_LHSIdIndices.GetSize(), ' '), rhs_tmp( m_RHSIdIndices.GetSize(), ' ');

      for( size_t row( 0), n_rows( lhs_dataset_ptr->GetIdsReference().GetNumberRows()); row < n_rows; ++row)
      {
        // copy the desired character into the temp string
        const char *row_ptr( lhs_dataset_ptr->GetIdsReference()[ row]);
        for( size_t col( 0); col < lhs_n_id_cols; ++col)
        {
          lhs_tmp[ col] = row_ptr[ m_LHSIdIndices( col)];
        }
        lhs_map_property_to_index[ util::TrimString( lhs_tmp)] = storage::Pair< std::string, size_t>( lhs_tmp, row);
      }
      for( size_t row( 0), n_rows( rhs_dataset_ptr->GetIdsReference().GetNumberRows()); row < n_rows; ++row)
      {
        // copy the desired character into the temp string
        const char *row_ptr( rhs_dataset_ptr->GetIdsReference()[ row]);
        for( size_t col( 0); col < rhs_n_id_cols; ++col)
        {
          rhs_tmp[ col] = row_ptr[ m_RHSIdIndices( col)];
        }
        rhs_map_property_to_index[ util::TrimString( rhs_tmp)] = storage::Pair< std::string, size_t>( rhs_tmp, row);
      }

      // initialize a new data set
      util::ShPtr< descriptor::Dataset> dataset
      (
        new descriptor::Dataset
        (
          id_pairs_to_result.GetSize(),
          lhs_dataset_ptr->GetFeatureSize() + rhs_dataset_ptr->GetFeatureSize(),
          result_size,
          m_RHSIdIndices.GetSize() + m_LHSIdIndices.GetSize()
        )
      );

      for( size_t i( 0); i < id_pairs_to_result.GetSize(); ++i)
      {
        auto row( dataset->GetFeaturesReference().GetRow( i));
        auto &lhs_index_buffstring( lhs_map_property_to_index.GetValue( id_pairs_to_result( i).First()));
        row.CreateSubVectorReference( lhs_dataset_ptr->GetFeatureSize(), 0).CopyValues( lhs_dataset_ptr->GetFeaturesReference().GetRow( lhs_index_buffstring.Second()));
        auto &rhs_index_buffstring( rhs_map_property_to_index.GetValue( id_pairs_to_result( i).Second()));
        row.CreateSubVectorReference( rhs_dataset_ptr->GetFeatureSize(), lhs_dataset_ptr->GetFeatureSize()).CopyValues
        (
          rhs_dataset_ptr->GetFeaturesReference().GetRow( rhs_index_buffstring.Second())
        );
        auto result_row( dataset->GetResultsReference().GetRow( i));
        result_row.CopyValues( id_pairs_to_result( i).Third());
      }

      // return the generated data set
      return dataset;
    } // GenerateDataSet

    //! @brief get the size of the complete dataset without chunking
    //! @return the size of the complete dataset without chunking
    size_t RetrieveDataSetJoin::GetTotalDatasetSize() const
    {
      // read the data set from the file
      io::IFStream input;
      io::File::MustOpenIFStream( input, m_IdPairsFileName);

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
    size_t RetrieveDataSetJoin::GetNominalSize() const
    {
      return GetTotalDatasetSize();
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer RetrieveDataSetJoin::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription( "Join two datasets with a pair-specific result value");
      parameters.AddInitializer
      (
        "id pairs filename",
        "fixed-width label for each dataset pair (e.g. for joining protein-ligand datasets this may indicate the target)",
        io::Serialization::GetAgent( &m_IdPairsFileName)
      );

      parameters.AddInitializer
      (
        "lhs",
        "dataset retriever to call to get the LHS",
        io::Serialization::GetAgent( &m_LHS)
      );

      parameters.AddInitializer
      (
        "lhs id label",
        "LHS data label",
        io::Serialization::GetAgent( &m_LHSLabel)
      );

      parameters.AddInitializer
      (
        "rhs",
        "dataset retriever to call to get the RHS",
        io::Serialization::GetAgent( &m_RHS)
      );

      parameters.AddInitializer
      (
        "rhs id label",
        "RHS data label",
        io::Serialization::GetAgent( &m_RHSLabel)
      );

      parameters.AddInitializer
      (
        "result label",
        "Result data label",
        io::Serialization::GetAgent( &m_ResultLabel)
      );

      return parameters;
    } // GetParameters

    void RetrieveDataSetJoin::UpdateIDsIndex()
    {
      if( m_LHSIdIndices.IsEmpty())
      {
        m_LHSIdIndices = m_LHS->GetIdCodeWithSizes().GetPropertyIndices( m_LHSLabel);
      }
      if( m_RHSIdIndices.IsEmpty())
      {
        m_RHSIdIndices = m_RHS->GetIdCodeWithSizes().GetPropertyIndices( m_RHSLabel);
      }

    }

    bool RetrieveDataSetJoin::ReadInitializerSuccessHook
    (
      const util::ObjectDataLabel &LABEL,
      std::ostream &ERR_STREAM
    )
    {
      return true;
    }

  } // namespace model
} // namespace bcl
