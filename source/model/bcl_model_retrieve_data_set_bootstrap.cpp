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
#include "model/bcl_model_retrieve_data_set_bootstrap.h"

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
    const util::SiPtr< const util::ObjectInterface> RetrieveDataSetBootstrap::s_Instance
    (
      util::Enumerated< RetrieveDataSetBase>::AddInstance( new RetrieveDataSetBootstrap())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief Clone is the virtual Clone constructor
    //! @return a pointer to new LocatorAtomFurthest which is a copy of this
    RetrieveDataSetBootstrap *RetrieveDataSetBootstrap::Clone() const
    {
      return new RetrieveDataSetBootstrap( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &RetrieveDataSetBootstrap::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &RetrieveDataSetBootstrap::GetAlias() const
    {
      static const std::string s_Name( "Bootstrap");
      return s_Name;
    }

    //! @brief Set the code / label for the feature (1st part) of the data set
    //! @param CODE the new code
    //! @return the code / label for the feature (1st part) of the data set
    void RetrieveDataSetBootstrap::SelectFeatures( const util::ObjectDataLabel &CODE)
    {
      RetrieveDataSetBase::SelectFeatures( CODE);
      m_Implementation->SelectFeatures( CODE);
    }

    //! @brief Set the code / label for the result (2nd part) of the data set
    //! @return the code / label for the result (2nd part) of the data set
    void RetrieveDataSetBootstrap::SelectResults( const util::ObjectDataLabel &CODE)
    {
      RetrieveDataSetBase::SelectResults( CODE);
      m_Implementation->SelectResults( CODE);
    }

    //! @brief Set which id columns to retrieve
    //! @param CODE the id column names to retrieve
    void RetrieveDataSetBootstrap::SelectIds( const util::ObjectDataLabel &CODE)
    {
      RetrieveDataSetBase::SelectIds( CODE);
      m_Implementation->SelectIds( CODE);
    }

    //! @brief Get the code / label set for the feature (1st part) of the data set with sizes of each property
    //! @return the code / label for the feature (1st part) of the data set with sizes of each property
    //! the feature code set
    FeatureLabelSet RetrieveDataSetBootstrap::GetFeatureLabelsWithSizes() const
    {
      return m_Implementation->GetFeatureLabelsWithSizes();
    }

    //! @brief Get the code / label for the result (2nd part) of the data set with sizes of each property
    //! @return the code / label for the result (2nd part) of the data set with sizes of each property
    //! the feature code set
    FeatureLabelSet RetrieveDataSetBootstrap::GetResultCodeWithSizes() const
    {
      return m_Implementation->GetResultCodeWithSizes();
    }

    //! @brief Get the code / label for the ids of the data set with sizes of each property
    //! @return the code / label for the ids of the data set with sizes of each property
    //! the feature code set
    FeatureLabelSet RetrieveDataSetBootstrap::GetIdCodeWithSizes() const
    {
      return m_Implementation->GetIdCodeWithSizes();
    }

    //! @brief get the number of partitions requested by the user, along with the partition ids
    //! @return the number of partitions requested by the user, along with the partition ids
    storage::Pair< size_t, math::RangeSet< size_t> > RetrieveDataSetBootstrap::GetNumberPartitionsAndIds() const
    {
      return m_Implementation->GetNumberPartitionsAndIds();
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief generate dataset
    //! @return generated dataset
    util::ShPtr< descriptor::Dataset>
    RetrieveDataSetBootstrap::GenerateDataSet()
    {
      util::ShPtr< descriptor::Dataset> dataset( m_Implementation->GenerateDataSet());

      random::UniformDistribution uni;
      uni.SetSeed( random::GetGlobalRandom().GetFlagRandomSeed()->GetFirstParameter()->GetNumericalValue< size_t>());
      const size_t mx( dataset->GetSize() - 1);
      util::ShPtr< descriptor::Dataset> dataset_new( dataset->HardCopy());
      for( size_t src_row_id( 0); src_row_id <= mx; ++src_row_id)
      {
        const size_t row_id( uni.Random( size_t( 0), mx));
        dataset_new->GetFeaturesReference().GetRow( src_row_id).CopyValues( dataset->GetFeaturesReference().GetRow( row_id));
        dataset_new->GetResultsReference().GetRow( src_row_id).CopyValues( dataset->GetResultsReference().GetRow( row_id));
        dataset_new->GetIdsReference().GetRow( src_row_id).CopyValues( dataset->GetIdsReference().GetRow( row_id));
      }
      return dataset_new;
    } // GenerateDataSet

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer RetrieveDataSetBootstrap::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription( "Retrieves a bootstrapped dataset. This means that the data are loaded in,");
      parameters.AddInitializer
      (
        "",
        "dataset retriever to call to get the entire data set",
        io::Serialization::GetAgent( &m_Implementation)
      );

      return parameters;
    } // GetParameters

    //! @brief get the nominal (e.g. best estimate without generating the entire dataset) size of the dataset
    //! @return the nominal (e.g. best estimate without generating the entire dataset) size of the dataset
    size_t RetrieveDataSetBootstrap::GetNominalSize() const
    {
      return GetNominalSize();
    }

  } // namespace model
} // namespace bcl
