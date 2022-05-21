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
#include "model/bcl_model_retrieve_data_set_balanced.h"

// includes from bcl - sorted alphabetically
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
    const util::SiPtr< const util::ObjectInterface> RetrieveDataSetBalanced::s_Instance
    (
      util::Enumerated< RetrieveDataSetBase>::AddInstance( new RetrieveDataSetBalanced())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief Clone is the virtual Clone constructor
    //! @return a pointer to new LocatorAtomFurthest which is a copy of this
    RetrieveDataSetBalanced *RetrieveDataSetBalanced::Clone() const
    {
      return new RetrieveDataSetBalanced( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &RetrieveDataSetBalanced::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief Set the code / label for the feature (1st part) of the data set
    //! @param CODE the new code
    //! @return the code / label for the feature (1st part) of the data set
    void RetrieveDataSetBalanced::SelectFeatures( const util::ObjectDataLabel &CODE)
    {
      RetrieveDataSetBase::SelectFeatures( CODE);
      for
      (
        storage::Vector< util::Implementation< RetrieveDataSetBase> >::iterator
          itr( m_Retrievers.Begin()), itr_end( m_Retrievers.End());
        itr != itr_end;
        ++itr
      )
      {
        // set the feature code for all the implementations
        ( *itr)->SelectFeatures( CODE);
      }
    }

    //! @brief Set the code / label for the result (2nd part) of the data set
    //! @return the code / label for the result (2nd part) of the data set
    void RetrieveDataSetBalanced::SelectResults( const util::ObjectDataLabel &CODE)
    {
      RetrieveDataSetBase::SelectResults( CODE);
      for
      (
        storage::Vector< util::Implementation< RetrieveDataSetBase> >::iterator
          itr( m_Retrievers.Begin()), itr_end( m_Retrievers.End());
        itr != itr_end;
        ++itr
      )
      {
        // set the result code for all the implementations
        ( *itr)->SelectResults( CODE);
      }
    }

    //! @brief Set which id columns to retrieve
    //! @param CODE the id column names to retrieve
    void RetrieveDataSetBalanced::SelectIds( const util::ObjectDataLabel &CODE)
    {
      RetrieveDataSetBase::SelectIds( CODE);
      for
      (
        storage::Vector< util::Implementation< RetrieveDataSetBase> >::iterator
          itr( m_Retrievers.Begin()), itr_end( m_Retrievers.End());
        itr != itr_end;
        ++itr
      )
      {
        // set the id code for all the implementations
        ( *itr)->SelectIds( CODE);
      }
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &RetrieveDataSetBalanced::GetAlias() const
    {
      static const std::string s_Name( "Balanced");
      return s_Name;
    }

    //! @brief Get the code / label set for the feature (1st part) of the data set with sizes of each property
    //! @return the code / label for the feature (1st part) of the data set with sizes of each property
    //! the feature code set
    FeatureLabelSet RetrieveDataSetBalanced::GetFeatureLabelsWithSizes() const
    {
      // test for null dataset
      if( m_Retrievers.IsEmpty())
      {
        return FeatureLabelSet();
      }

      // otherwise, get the feature label set from the first retriever (they should both be the same)
      return m_Retrievers.FirstElement()->GetFeatureLabelsWithSizes();
    }

    //! @brief Get the code / label for the result (2nd part) of the data set with sizes of each property
    //! @return the code / label for the result (2nd part) of the data set with sizes of each property
    //! the feature code set
    FeatureLabelSet RetrieveDataSetBalanced::GetResultCodeWithSizes() const
    {
      // test for null dataset
      if( m_Retrievers.IsEmpty())
      {
        return FeatureLabelSet();
      }

      // otherwise, get the feature label set from the first retriever (they should both be the same)
      return m_Retrievers.FirstElement()->GetResultCodeWithSizes();
    }

    //! @brief Get the code / label for the ids of the data set with sizes of each property
    //! @return the code / label for the ids of the data set with sizes of each property
    //! the feature code set
    FeatureLabelSet RetrieveDataSetBalanced::GetIdCodeWithSizes() const
    {
      // test for null dataset
      if( m_Retrievers.IsEmpty())
      {
        return FeatureLabelSet();
      }

      // otherwise, get the feature label set from the first retriever (they should both be the same)
      return m_Retrievers.FirstElement()->GetIdCodeWithSizes();
    }

    //! @brief get whether dataset generation requires labels
    //! @return true if dataset generation requires labels
    bool RetrieveDataSetBalanced::RequiresFeatureLabels() const
    {
      // iterate through the retrievers
      for
      (
        storage::Vector< util::Implementation< RetrieveDataSetBase> >::const_iterator
          itr( m_Retrievers.Begin()), itr_end( m_Retrievers.End());
        itr != itr_end;
        ++itr
      )
      {
        // determine whether this retriever requires feature labels
        if( ( *itr)->RequiresFeatureLabels())
        {
          // one retriever requires the feature labels, so return true
          return true;
        }
      }

      return false;
    }

    //! @brief get whether dataset generation requires result labels
    //! @return true if dataset generation requires result labels
    bool RetrieveDataSetBalanced::RequiresResultLabels() const
    {
      // iterate through the retrievers
      for
      (
        storage::Vector< util::Implementation< RetrieveDataSetBase> >::const_iterator
          itr( m_Retrievers.Begin()), itr_end( m_Retrievers.End());
        itr != itr_end;
        ++itr
      )
      {
        // determine whether this retriever requires result labels
        if( ( *itr)->RequiresResultLabels())
        {
          // one retriever requires the result labels, so return true
          return true;
        }
      }

      return false;
    }

    //! @brief get the number of partitions requested by the user, along with the partition ids
    //! @return the number of partitions requested by the user, along with the partition ids
    storage::Pair< size_t, math::RangeSet< size_t> > RetrieveDataSetBalanced::GetNumberPartitionsAndIds() const
    {
      if( m_Retrievers.IsEmpty())
      {
        return
          storage::Pair< size_t, math::RangeSet< size_t> >( 1, math::RangeSet< size_t>( math::Range< size_t>( 0, 0)));
      }
      storage::Pair< size_t, math::RangeSet< size_t> > initial
      (
        m_Retrievers.FirstElement()->GetNumberPartitionsAndIds()
      );
      for
      (
        storage::Vector< util::Implementation< RetrieveDataSetBase> >::const_iterator
          itr( m_Retrievers.Begin() + 1), itr_end( m_Retrievers.End());
        itr != itr_end;
        ++itr
      )
      {
        // determine whether this retriever requires result labels
        if( ( *itr)->GetNumberPartitionsAndIds() != initial)
        {
          // different partitions, return undefined
          return
            storage::Pair< size_t, math::RangeSet< size_t> >( util::GetUndefined< size_t>(), math::RangeSet< size_t>());
        }
      }
      // all partitions are the same, return it
      return initial;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief generate dataset
    //! @return generated dataset
    util::ShPtr< descriptor::Dataset>
    RetrieveDataSetBalanced::GenerateDataSet()
    {
      // just return if there were no datasets
      if( m_Retrievers.IsEmpty())
      {
        return util::ShPtr< descriptor::Dataset>( new descriptor::Dataset());
      }
      else if( m_Retrievers.GetSize() == 1) // if there was only 1 retriever, just return its results
      {
        return m_Retrievers.FirstElement()->GenerateDataSet();
      }

      // initialize a new vector of data sets, beginning with the first element
      util::ShPtrVector< descriptor::Dataset> datasets
      (
        1, m_Retrievers.FirstElement()->GenerateDataSet()
      );

      // determine feature / result size to ensure that all features and results are of the same size
      const size_t feature_size( datasets.FirstElement()->GetFeatureSize());
      const size_t result_size( datasets.FirstElement()->GetResultSize());
      const size_t id_size( datasets.FirstElement()->GetIdSize());

      // keep track of the size of the largest data set
      size_t largest_dataset_size( datasets.FirstElement()->GetSize());

      // generate the rest of the datasets
      for
      (
        size_t dataset_count( 1), number_datasets( m_Retrievers.GetSize());
        dataset_count < number_datasets;
        ++dataset_count
      )
      {
        // create a data set using the next retriever
        datasets.PushBack( m_Retrievers( dataset_count)->GenerateDataSet());

        // ensure that the feature size and result size were correct
        BCL_Assert
        (
          feature_size == datasets.LastElement()->GetFeatureSize(),
          " cannot balance datasets with differing feature sizes"
        );
        BCL_Assert
        (
          result_size == datasets.LastElement()->GetResultSize(),
          " cannot balance datasets with differing result sizes"
        );
        BCL_Assert
        (
          id_size == datasets.LastElement()->GetIdSize(),
          " cannot balance datasets with differing id sizes"
        );

        // update the largest dataset size
        largest_dataset_size = std::max( largest_dataset_size, datasets.LastElement()->GetSize());
      }

      const size_t dataset_size( largest_dataset_size * datasets.GetSize());

      // create a new data set memory that balances the output of the old data sets
      // initialize a new data set
      util::ShPtr< descriptor::Dataset> sp_dataset_balanced
      (
        new descriptor::Dataset
        (
          dataset_size,
          GetFeatureLabelsWithSizes(),
          GetResultCodeWithSizes(),
          GetIdCodeWithSizes()
        )
      );
      storage::Vector< linal::MatrixConstReference< float> > unbalanced_features( datasets.GetSize());
      storage::Vector< linal::MatrixConstReference< float> > unbalanced_results( datasets.GetSize());
      storage::Vector< linal::MatrixConstReference< char> > unbalanced_ids( datasets.GetSize());
      for
      (
        size_t dataset_number( 0), number_datasets( datasets.GetSize());
        dataset_number < number_datasets;
        ++dataset_number
      )
      {
        unbalanced_features( dataset_number) = linal::MatrixConstReference< float>( datasets( dataset_number)->GetFeaturesReference());
        unbalanced_results( dataset_number)  = linal::MatrixConstReference< float>( datasets( dataset_number)->GetResultsReference() );
        unbalanced_ids( dataset_number)      = linal::MatrixConstReference< char>( datasets( dataset_number)->GetIdsReference()     );
      }
      for( size_t data_number( 0), balanced_row_number( 0); data_number < largest_dataset_size; ++data_number)
      {
        for
        (
          size_t dataset_number( 0), number_datasets( datasets.GetSize());
          dataset_number < number_datasets;
          ++dataset_number, ++balanced_row_number
        )
        {
          // get a reference to the dataset that this feature/result will be pulled from
          const descriptor::Dataset &dataset( *datasets( dataset_number));

          // determine the feature number to take from this dataset
          const size_t feature_number( data_number % dataset.GetSize());

          sp_dataset_balanced->AddData
          (
            balanced_row_number,
            unbalanced_features( dataset_number).GetRow( feature_number),
            unbalanced_results( dataset_number).GetRow( feature_number),
            unbalanced_ids( dataset_number).GetRow( feature_number)
          );
        }
      }

      return sp_dataset_balanced;
    } // GenerateDataSet

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer RetrieveDataSetBalanced::GetSerializer() const
    {
      io::Serializer member_data;
      member_data.SetClassDescription
      (
        "balances features/results from several data sets such that a randomly chosen element in the data set "
        "has a roughly equal probability of coming from each data set"
      );
      member_data.AddInitializer( "", "datasets to balance", io::Serialization::GetAgent( &m_Retrievers));

      return member_data;
    } // GetParameters

    //! @brief get the nominal (e.g. best estimate without generating the entire dataset) size of the dataset
    //! @return the nominal (e.g. best estimate without generating the entire dataset) size of the dataset
    size_t RetrieveDataSetBalanced::GetNominalSize() const
    {
      size_t max_size( 0);
      for
      (
        storage::Vector< util::Implementation< RetrieveDataSetBase> >::const_iterator
          itr( m_Retrievers.Begin()), itr_end( m_Retrievers.End());
        itr != itr_end;
        ++itr
      )
      {
        // set the result code for all the implementations
        max_size = std::max( max_size, ( *itr)->GetNominalSize());
      }
      return max_size * m_Retrievers.GetSize();
    }

  } // namespace model
} // namespace bcl
