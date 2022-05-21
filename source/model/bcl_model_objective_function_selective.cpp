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
#include "model/bcl_model_objective_function_selective.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"
#include "math/bcl_math_running_average.h"
#include "model/bcl_model_data_set_select_columns.h"
#include "model/bcl_model_feature_label_set.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {

    // add the interface to the set of known implementations
    const util::SiPtr< const util::ObjectInterface> ObjectiveFunctionSelective::s_Instance
    (
      util::Enumerated< ObjectiveFunctionInterface>::AddInstance
      (
        new ObjectiveFunctionSelective()
      )
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! copy constructor
    ObjectiveFunctionSelective *ObjectiveFunctionSelective::Clone() const
    {
      return new ObjectiveFunctionSelective( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &ObjectiveFunctionSelective::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief determine what sign of the derivative of this objective function indicates improvement
    //! @return the sign of the derivative of this objective function indicates improvement
    opti::ImprovementType ObjectiveFunctionSelective::GetImprovementType() const
    {
      if( m_Objective.IsDefined())
      {
        return m_Objective->GetImprovementType();
      }
      return opti::e_LargerIsBetter;
    }

    //! @brief get the overall goal of the objective function
    //! @return the goal of the objective function
    ObjectiveFunctionInterface::Goal ObjectiveFunctionSelective::GetGoalType() const
    {
      if( m_Objective.IsDefined())
      {
        return m_Objective->GetGoalType();
      }
      return e_Other;
    }

    //! @brief classify. Obtain a matrix of with N|n for all predicted negatives, P|p for all predicted positives, \0 for all non-predicted values
    //!        case indicates whether the prediction was true (upper) or false (lower)
    //! @param EXPERIMENTAL feature dataset with experimental values
    //! @param PREDICTED feature dataset with predicted values
    //! @return matrix with PpNn\0 values indicating TP,FP,TN,FN,NA
    linal::Matrix< char> ObjectiveFunctionSelective::GetFeaturePredictionClassifications
    (
      const FeatureDataSetInterface< float> &EXPERIMENTAL,
      const FeatureDataSetInterface< float> &PREDICTED
    ) const
    {
      return m_Objective->GetFeaturePredictionClassifications( EXPERIMENTAL, PREDICTED);
    }

    //! @brief set dataset; some objective functions may need to setup internal data structures based on this function
    //! @param DATA monitoring dataset results, non-scaled
    //! @param IDS ids; can be used by the objective function
    void ObjectiveFunctionSelective::SetData
    (
      const FeatureDataSet< float> &DATA,
      const FeatureDataSet< char> &IDS
    )
    {
      m_MappingIds.Reset();

      // get the feature labels so that sequences can be identified
      DataSetSelectColumns sequence_id_selector;
      if( !m_DbIdLabel.IsEmpty())
      {
        const FeatureLabelSet model_feature_labels( *IDS.GetFeatureLabelSet());
        sequence_id_selector =
          DataSetSelectColumns
          (
            model_feature_labels.GetSize(),
            model_feature_labels.GetPropertyIndices
            (
              util::ObjectDataLabel( m_DbIdLabel.GetValue(), m_DbIdLabel.GetArguments())
            )
          );
      }
      else
      {
        // no labels were specified, use them all
        sequence_id_selector =
          DataSetSelectColumns( IDS.GetFeatureSize(), storage::CreateIndexVector( IDS.GetFeatureSize()));
      }

      BCL_Assert( sequence_id_selector.GetOutputFeatureSize(), "Could not find " + m_DbIdLabel.ToString() + " in IDs!");

      const size_t dataset_size( DATA.GetNumberFeatures());
      BCL_Assert( dataset_size == IDS.GetNumberFeatures(), "IDs and Data have different size!");

      {
        std::string this_sequence_id( sequence_id_selector.GetOutputFeatureSize(), ' ');
        for( size_t i( 0); i < dataset_size; ++i)
        {
          sequence_id_selector( IDS( i), &this_sequence_id[ 0]);
          m_MappingIds[ this_sequence_id].PushBack( i);
        }
      }

      m_ExperimentalPruned = FeatureDataSet< float>( m_MappingIds.GetSize(), DATA.GetFeatureSize(), float( 0.0));
      FeatureDataSet< char> pruned_ids( m_MappingIds.GetSize(), IDS.GetFeatureSize(), ' ');
      size_t molecule_index( 0);
      for
      (
        storage::Map< std::string, storage::Vector< size_t> >::const_iterator
          itr( m_MappingIds.Begin()), itr_end( m_MappingIds.End());
          itr != itr_end;
        ++itr, ++molecule_index
      )
      {
        m_ExperimentalPruned.GetRawMatrix().ReplaceRow( molecule_index, DATA( itr->second.FirstElement()));
        pruned_ids.GetRawMatrix().ReplaceRow( molecule_index, IDS( itr->second.FirstElement()));
      }
      m_Objective->SetData( m_ExperimentalPruned, pruned_ids);
    }

    //! @brief evaluate for a given datasets with experimental and predicted values the implemented objective function
    //! @param EXPERIMENTAL feature dataset with experimental values
    //! @param PREDICTED feature dataset with predicted values
    //! @return objective function value based on the given data of experimental and predicted values
    float ObjectiveFunctionSelective::operator()
    (
      const FeatureDataSetInterface< float> &EXPERIMENTAL,
      const FeatureDataSetInterface< float> &PREDICTED
    ) const
    {
      FeatureDataSet< float> predicted( m_MappingIds.GetSize(), EXPERIMENTAL.GetFeatureSize(), float( 0.0));

      size_t molecule_index( 0);
      math::RunningAverage< linal::Vector< float> > minmax;
      for
      (
        storage::Map< std::string, storage::Vector< size_t> >::const_iterator
          itr( m_MappingIds.Begin()), itr_end( m_MappingIds.End());
          itr != itr_end;
        ++itr, ++molecule_index
      )
      {
        minmax.Reset();
        const storage::Vector< size_t> &values( itr->second);
        for
        (
          storage::Vector< size_t>::const_iterator itr_vec( values.Begin()), itr_vec_end( values.End());
          itr_vec != itr_vec_end;
          ++itr_vec
        )
        {
          minmax += PREDICTED( *itr_vec);
        }
        predicted.GetRawMatrix().ReplaceRow( molecule_index, minmax.GetAverage());
      }

      return m_Weight * m_Objective->operator ()( m_ExperimentalPruned, predicted);
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer ObjectiveFunctionSelective::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Filters results; selecting only the maximum value among results with the same ID"
      );

      parameters.AddInitializer
      (
        "weight",
        "amount by which to weight the objective function's output",
        io::Serialization::GetAgent( &m_Weight),
        "1.0"
      );

      parameters.AddInitializer
      (
        "function",
        "core objective function to use",
        io::Serialization::GetAgent( &m_Objective)
      );

      parameters.AddInitializer
      (
        "label",
        "id label. If not specified, uses the full id label from the dataset",
        io::Serialization::GetAgent( &m_DbIdLabel),
        ""
      );

      return parameters;
    }

  } // namespace model
} // namespace bcl
