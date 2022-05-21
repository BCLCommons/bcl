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
#include "model/bcl_model_retrieve_data_set_combined.h"

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
    const util::SiPtr< const util::ObjectInterface> RetrieveDataSetCombined::s_Instance
    (
      util::Enumerated< RetrieveDataSetBase>::AddInstance( new RetrieveDataSetCombined())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief Clone is the virtual Clone constructor
    //! @return a pointer to new LocatorAtomFurthest which is a copy of this
    RetrieveDataSetCombined *RetrieveDataSetCombined::Clone() const
    {
      return new RetrieveDataSetCombined( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &RetrieveDataSetCombined::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief Set the code / label for the feature (1st part) of the data set
    //! @param CODE the new code
    //! @return the code / label for the feature (1st part) of the data set
    void RetrieveDataSetCombined::SelectFeatures( const util::ObjectDataLabel &CODE)
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
    void RetrieveDataSetCombined::SelectResults( const util::ObjectDataLabel &CODE)
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
    void RetrieveDataSetCombined::SelectIds( const util::ObjectDataLabel &CODE)
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
    const std::string &RetrieveDataSetCombined::GetAlias() const
    {
      static const std::string s_Name( "Combined");
      return s_Name;
    }

    //! @brief Get the code / label set for the feature (1st part) of the data set with sizes of each property
    //! @return the code / label for the feature (1st part) of the data set with sizes of each property
    //! the feature code set
    FeatureLabelSet RetrieveDataSetCombined::GetFeatureLabelsWithSizes() const
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
    FeatureLabelSet RetrieveDataSetCombined::GetResultCodeWithSizes() const
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
    FeatureLabelSet RetrieveDataSetCombined::GetIdCodeWithSizes() const
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
    bool RetrieveDataSetCombined::RequiresFeatureLabels() const
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
    bool RetrieveDataSetCombined::RequiresResultLabels() const
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

  ////////////////
  // operations //
  ////////////////

    //! @brief generate dataset
    //! @return generated dataset
    util::ShPtr< descriptor::Dataset>
    RetrieveDataSetCombined::GenerateDataSet()
    {
      // determine the approximate size of the complete dataset
      const size_t full_dataset_size( GetNominalSize());

      // initialize a new data set
      util::ShPtr< descriptor::Dataset> dataset
      (
        new descriptor::Dataset
        (
          full_dataset_size,
          GetFeatureLabelsWithSizes(),
          GetResultCodeWithSizes(),
          GetIdCodeWithSizes()
        )
      );

      linal::MatrixReference< float> features( dataset->GetFeaturesReference());
      linal::MatrixReference< float> results( dataset->GetResultsReference());
      linal::MatrixReference< char> ids( dataset->GetIdsReference());

      size_t number_features_so_far( 0);
      for
      (
        storage::Vector< util::Implementation< RetrieveDataSetBase> >::iterator
          itr( m_Retrievers.Begin()), itr_end( m_Retrievers.End());
        itr != itr_end;
        ++itr
      )
      {
        // get the nominal size of this dataset
        const size_t nominal_dataset_size( ( *itr)->GetNominalSize());

        // put the data set from the next retriever directly into the dataset's matrices
        number_features_so_far +=
          ( *itr)->GenerateDataSubset
          (
            math::Range< size_t>
            (
              math::RangeBorders::e_LeftClosed,
              0,
              nominal_dataset_size,
              math::RangeBorders::e_RightOpen
            ),
            features,
            results,
            ids,
            number_features_so_far
          );
      }

      // remove unused rows
      dataset->ShrinkRows( number_features_so_far);

      // return the generated data set
      return dataset;
    } // GenerateDataSet

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer RetrieveDataSetCombined::GetSerializer() const
    {
      io::Serializer member_data;
      member_data.SetClassDescription( "concatenate data sets retrieved from each implementation");
      member_data.AddInitializer( "", "datasets to combine", io::Serialization::GetAgent( &m_Retrievers));

      return member_data;
    } // GetParameters

    //! @brief test whether this retriever can generate sub-ranges of datasets without loading the entire dataset
    //! @return true if this retriever can generate sub-ranges of datasets without loading the entire dataset
    bool RetrieveDataSetCombined::SupportsEfficientSubsetLoading() const
    {
      // test all internal datasets

      // test that each data set retriever supports efficient subset loading
      for
      (
        storage::Vector< util::Implementation< RetrieveDataSetBase> >::const_iterator
          itr( m_Retrievers.Begin()), itr_end( m_Retrievers.End());
        itr != itr_end;
        ++itr
      )
      {
        if( !( *itr)->SupportsEfficientSubsetLoading())
        {
          return false;
        }
      }

      // all internal implementations support efficient subset loading and GetSize
      return true;
    }

    //! @brief get the number of partitions requested by the user, along with the partition ids
    //! @return the number of partitions requested by the user, along with the partition ids
    storage::Pair< size_t, math::RangeSet< size_t> > RetrieveDataSetCombined::GetNumberPartitionsAndIds() const
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

    //! @brief load a range of data from the dataset
    //! @param SUBSET the range of data to load
    //! @param FEATURES_STORAGE where to store features that are loaded, must be large enough to hold the subset without resizing
    //! @param RESULTS_STORAGE where to store the corresponding results, must be large enough to hold the subset without resizing
    //! @param START_FEATURE_NUMBER position to store the first feature in FEATURES_STORAGE
    //! @return # of features actually loaded
    //! Note: Implementations should overload this and SupportsEfficientSubsetLoading together
    size_t RetrieveDataSetCombined::GenerateDataSubset
    (
      const math::Range< size_t> &SUBSET,
      linal::MatrixInterface< float> &FEATURES_STORAGE,
      linal::MatrixInterface< float> &RESULTS_STORAGE,
      linal::MatrixInterface< char> &IDS_STORAGE,
      const size_t &START_FEATURE_NUMBER
    )
    {
      // ensure that features and results are the same size
      BCL_Assert
      (
        FEATURES_STORAGE.GetNumberRows() == RESULTS_STORAGE.GetNumberRows(),
        "Different number of features and results"
      );

      // ensure that the start position is within the storage matrix
      BCL_Assert
      (
        START_FEATURE_NUMBER <= FEATURES_STORAGE.GetNumberRows(),
        "start feature was higher than storage"
      );

      // standardize the range, so the left side of it has a closed border and the right has an open border
      const math::Range< size_t> standard_subset( SUBSET.StandardizeRange());

      // ensure that the end position is within the storage matrix
      BCL_Assert
      (
        START_FEATURE_NUMBER + standard_subset.GetWidth() <= FEATURES_STORAGE.GetNumberRows(),
        "Storage was too small to contain entire subset"
      );

      // standardize the range, so the left side of it has a closed border and the right has an open border
      size_t first_to_pick( standard_subset.GetMin());
      const size_t total_size( standard_subset.GetWidth());

      // keep track of the nominal and actual sizes
      // actual size may be smaller than nominal size, e.g. if some features could not be generated or were
      // of the wrong size
      size_t nominal_size_so_far( 0), actual_size_so_far( 0);

      // get the size from each of the internal datasets
      for
      (
        storage::Vector< util::Implementation< RetrieveDataSetBase> >::iterator
          itr( m_Retrievers.Begin()), itr_end( m_Retrievers.End());
        itr != itr_end && nominal_size_so_far < total_size;
        ++itr
      )
      {
        // use GetSize to decide whether to generate this dataset or not
        size_t size( ( *itr)->GetNominalSize());

        // if we have added any elements to the dataset yet, then continue adding them
        if( nominal_size_so_far != 0)
        {
          // just append the subset of interest
          actual_size_so_far +=
            ( *itr)->GenerateDataSubset
            (
              math::Range< size_t>
              (
                math::RangeBorders::e_LeftClosed,
                0,
                total_size - nominal_size_so_far,
                math::RangeBorders::e_RightOpen
              ),
              FEATURES_STORAGE,
              RESULTS_STORAGE,
              IDS_STORAGE,
              actual_size_so_far + START_FEATURE_NUMBER
            );
          nominal_size_so_far += size;
        }
        else
        {
          if( size <= first_to_pick)
          {
            // subtract size from first to pick and proceed to the next dataset
            first_to_pick -= size;
          }
          else
          {
            // store the new dataset directly in dataset
            actual_size_so_far =
              ( *itr)->GenerateDataSubset
              (
                math::Range< size_t>
                (
                  math::RangeBorders::e_LeftClosed,
                  first_to_pick,
                  std::min( size, total_size + first_to_pick),
                  math::RangeBorders::e_RightOpen
                ),
                FEATURES_STORAGE,
                RESULTS_STORAGE,
                IDS_STORAGE,
                START_FEATURE_NUMBER
              );
            nominal_size_so_far += size - first_to_pick;
            first_to_pick = 0;
          }
        }
      }

      return actual_size_so_far;
    }

    //! @brief get the nominal (e.g. best estimate without generating the entire dataset) size of the dataset
    //! @return the nominal (e.g. best estimate without generating the entire dataset) size of the dataset
    size_t RetrieveDataSetCombined::GetNominalSize() const
    {
      size_t size( 0);
      // get the size from each of the internal datasets
      for
      (
        storage::Vector< util::Implementation< RetrieveDataSetBase> >::const_iterator
          itr( m_Retrievers.Begin()), itr_end( m_Retrievers.End());
        itr != itr_end;
        ++itr
      )
      {
        size += ( *itr)->GetNominalSize();
      }

      return size;
    }

  } // namespace model
} // namespace bcl
