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

// include forward header of this class
#include "model/bcl_model_score_dataset_partition.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "model/bcl_model_approximator_decision_tree.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> ScoreDatasetPartition::s_Instance
    (
      util::Enumerated< ScoreDatasetInterface>::AddInstance( new ScoreDatasetPartition())
    );

    //! @brief Clone function
    //! @return pointer to new ScoreDatasetPartition
    ScoreDatasetPartition *ScoreDatasetPartition::Clone() const
    {
      return new ScoreDatasetPartition( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &ScoreDatasetPartition::GetClassIdentifier() const
    {
      return GetStaticClassName< ScoreDatasetPartition>();
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &ScoreDatasetPartition::GetAlias() const
    {
      static const std::string s_Name( "Partition");
      return s_Name;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief score a given dataset
    //! @param DATASET dataset of interest
    //! @return scores of the dataset
    linal::Vector< float> ScoreDatasetPartition::Score( const descriptor::Dataset &DATASET) const
    {
      // create an iterate to setup the feature result and state vector
      ApproximatorDecisionTree iterate
      (
        m_Partitioner,
        util::ShPtr< ObjectiveFunctionWrapper>(),
        m_Cutoff,
        util::ShPtr< descriptor::Dataset>()
      );

      // create the data set references
      util::ShPtr< storage::Vector< FeatureResultAndState> > processed_dataset
      (
        iterate.CreateDataSetReferences( DATASET.GetFeaturesReference(), DATASET.GetResultsReference())
      );

      return m_Partitioner->GetAllPartitionRatings( *processed_dataset, true);
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer ScoreDatasetPartition::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "uses the split rating from a decision tree dataset partitioner"
      );

      parameters.AddInitializer
      (
        "cutoff",
        "result value that separates actives from inactives",
        io::Serialization::GetAgent( &m_Cutoff),
        "0"
      );

      parameters.AddInitializer
      (
        "partitioner",
        "method of deciding which component of the feature vector to use to add a new node to the decision tree",
        io::Serialization::GetAgent( &m_Partitioner),
        "InformationGain"
      );

      return parameters;
    }

  } // namespace model
} // namespace bcl
