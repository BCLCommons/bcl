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
#include "model/bcl_model_retrieve_data_set_rescaled.h"

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
    const util::SiPtr< const util::ObjectInterface> RetrieveDataSetRescaled::s_Instance
    (
      util::Enumerated< RetrieveDataSetBase>::AddInstance( new RetrieveDataSetRescaled())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    RetrieveDataSetRescaled::RetrieveDataSetRescaled() :
      m_RescaleInputType( RescaleFeatureDataSet::e_None),
      m_RescaleInputRange( -1, 1),
      m_RescaleOutputType( RescaleFeatureDataSet::e_None),
      m_RescaleOutputRange( -1, 1)
    {
    }

    //! @brief Clone is the virtual Clone constructor
    //! @return a pointer to new LocatorAtomFurthest which is a copy of this
    RetrieveDataSetRescaled *RetrieveDataSetRescaled::Clone() const
    {
      return new RetrieveDataSetRescaled( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &RetrieveDataSetRescaled::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &RetrieveDataSetRescaled::GetAlias() const
    {
      static const std::string s_Name( "Rescale");
      return s_Name;
    }

    //! @brief Set the code / label for the feature (1st part) of the data set
    //! @param CODE the new code
    //! @return the code / label for the feature (1st part) of the data set
    void RetrieveDataSetRescaled::SelectFeatures( const util::ObjectDataLabel &CODE)
    {
      RetrieveDataSetBase::SelectFeatures( CODE);
      m_Implementation->SelectFeatures( CODE);
    }

    //! @brief Set the code / label for the result (2nd part) of the data set
    //! @return the code / label for the result (2nd part) of the data set
    void RetrieveDataSetRescaled::SelectResults( const util::ObjectDataLabel &CODE)
    {
      RetrieveDataSetBase::SelectResults( CODE);
      m_Implementation->SelectResults( CODE);
    }

    //! @brief Set which id columns to retrieve
    //! @param CODE the id column names to retrieve
    void RetrieveDataSetRescaled::SelectIds( const util::ObjectDataLabel &CODE)
    {
      RetrieveDataSetBase::SelectIds( CODE);
      m_Implementation->SelectIds( CODE);
    }

    //! @brief Get the code / label set for the feature (1st part) of the data set with sizes of each property
    //! @return the code / label for the feature (1st part) of the data set with sizes of each property
    //! the feature code set
    FeatureLabelSet RetrieveDataSetRescaled::GetFeatureLabelsWithSizes() const
    {
      return m_Implementation->GetFeatureLabelsWithSizes();
    }

    //! @brief Get the code / label for the result (2nd part) of the data set with sizes of each property
    //! @return the code / label for the result (2nd part) of the data set with sizes of each property
    //! the feature code set
    FeatureLabelSet RetrieveDataSetRescaled::GetResultCodeWithSizes() const
    {
      return m_Implementation->GetResultCodeWithSizes();
    }

    //! @brief Get the code / label for the ids of the data set with sizes of each property
    //! @return the code / label for the ids of the data set with sizes of each property
    //! the feature code set
    FeatureLabelSet RetrieveDataSetRescaled::GetIdCodeWithSizes() const
    {
      return m_Implementation->GetIdCodeWithSizes();
    }

    //! @brief get the number of partitions requested by the user, along with the partition ids
    //! @return the number of partitions requested by the user, along with the partition ids
    storage::Pair< size_t, math::RangeSet< size_t> > RetrieveDataSetRescaled::GetNumberPartitionsAndIds() const
    {
      return m_Implementation->GetNumberPartitionsAndIds();
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief generate dataset
    //! @return generated dataset
    util::ShPtr< descriptor::Dataset>
    RetrieveDataSetRescaled::GenerateDataSet()
    {
      util::ShPtr< descriptor::Dataset> dataset_sp( m_Implementation->GenerateDataSet());
      dataset_sp->GetFeatures().Rescale( m_RescaleInputRange, m_RescaleInputType);
      dataset_sp->GetResults().Rescale( m_RescaleOutputRange, m_RescaleOutputType);

      return dataset_sp;
    } // GenerateDataSet

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer RetrieveDataSetRescaled::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Rescales a dataset. This is *only* useful for combining datasets that have different logical units or scales. "
        "For training machine learning models, set the rescaling within the model itself"
      );
      parameters.AddInitializer
      (
        "",
        "dataset retriever to call to get the entire data set",
        io::Serialization::GetAgent( &m_Implementation)
      );
      parameters.AddInitializer
      (
        "input scaling",
        "Type of input scaling. Normally AveStd works best, but MinMax and None may also be used in some circumstances",
        io::Serialization::GetAgent( &m_RescaleInputType),
        "AveStd"
      );
      parameters.AddInitializer
      (
        "output scaling",
        "Type of input scaling. Normally AveStd works best, but MinMax and None may also be used in some circumstances",
        io::Serialization::GetAgent( &m_RescaleOutputType),
        "None"
      );
      parameters.AddInitializer
      (
        "rescale input to",
        "rescale range. All input features will be rescaled to this range unless input scaling=None",
        io::Serialization::GetAgent( &m_RescaleInputRange),
        "\"[-1,1]\""
      );
      parameters.AddInitializer
      (
        "rescale output to",
        "rescale range. All results will be rescaled to this range unless output scaling=None",
        io::Serialization::GetAgent( &m_RescaleOutputRange),
        "\"[-1,1]\""
      );

      return parameters;
    } // GetParameters

    //! @brief get the nominal (e.g. best estimate without generating the entire dataset) size of the dataset
    //! @return the nominal (e.g. best estimate without generating the entire dataset) size of the dataset
    size_t RetrieveDataSetRescaled::GetNominalSize() const
    {
      return m_Implementation->GetNominalSize();
    }

  } // namespace model
} // namespace bcl
