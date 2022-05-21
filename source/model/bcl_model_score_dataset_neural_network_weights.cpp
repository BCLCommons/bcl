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
#include "model/bcl_model_score_dataset_neural_network_weights.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"
#include "model/bcl_model_neural_network.h"
#include "util/bcl_util_sh_ptr_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> ScoreDatasetNeuralNetworkWeights::s_Instance
    (
      util::Enumerated< ScoreDatasetInterface>::AddInstance( new ScoreDatasetNeuralNetworkWeights())
    );

    // Map from initializer string/type to the shptrvector of model interfaces; saves loading potentially
    // gigantic models repetitively
    storage::Map< std::string, storage::Vector< util::OwnPtr< NeuralNetwork> > > ScoreDatasetNeuralNetworkWeights::s_Models =
      storage::Map< std::string, storage::Vector< util::OwnPtr< NeuralNetwork> > >();

    //! @brief Clone function
    //! @return pointer to new ScoreDatasetNeuralNetworkWeights
    ScoreDatasetNeuralNetworkWeights *ScoreDatasetNeuralNetworkWeights::Clone() const
    {
      return new ScoreDatasetNeuralNetworkWeights( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &ScoreDatasetNeuralNetworkWeights::GetClassIdentifier() const
    {
      return GetStaticClassName< ScoreDatasetNeuralNetworkWeights>();
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &ScoreDatasetNeuralNetworkWeights::GetAlias() const
    {
      static const std::string s_Name( "NeuralNetworkWeights");
      return s_Name;
    }

  ////////////////
  // operations //
  ////////////////

  ///////////////
  // operators //
  ///////////////

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief score features and results
    //! @param FEATURES matrix of features
    //! @param RESULTS matrix of results
    //! @return scores of the dataset
    linal::Vector< float> ScoreDatasetNeuralNetworkWeights::Score( const descriptor::Dataset &DATASET) const
    {
      const size_t feature_size( DATASET.GetFeatureSize());

      // handle empty dataset
      if( DATASET.GetSize() == size_t( 0))
      {
        return linal::Vector< float>( feature_size, float( 0.0));
      }

      const util::ShPtrVector< NeuralNetwork> &models( *m_Models);
      const size_t n_models( models.GetSize());

      // ensure that the # of features and results was correct
      BCL_Assert
      (
        models.FirstElement()->GetNumberInputs() == feature_size,
        "Please specify the same features as was used to train the networks"
      );

      // store all weight effect matrices
      storage::Vector< linal::Matrix< float> > weight_effects( n_models);

      // for each model
      // 1. Obtain its matrix of weights
      for( size_t model_number( 0); model_number < n_models; ++model_number)
      {
        // get the weights from this network
        const storage::Vector< linal::Matrix< float> > &weights( models( model_number)->GetWeight());
        // get the first layer weights
        linal::Matrix< float> weights_mult( weights( 0).Transposed());
        for( size_t layer_number( 1), n_layers( weights.GetSize()); layer_number < n_layers; ++layer_number)
        {
          weights_mult = weights_mult * weights( layer_number).Transposed();
        }
        weight_effects( model_number) = weights_mult;
      }

      // load information about the cross validation and how it performed
      storage::List< CrossValidationInfo> cv_info( m_ModelStorage->RetrieveEnsembleCVInfo());

      // find the set of models that belonged to each cross-validation independent set
      storage::Map< util::ObjectDataLabel, storage::Vector< float> > independent_to_result;
      storage::Map< util::ObjectDataLabel, float> independent_to_median;

      // iterate through the list of cv info
      for
      (
        storage::List< CrossValidationInfo>::const_iterator itr( cv_info.Begin()), itr_end( cv_info.End());
        itr != itr_end;
        ++itr
      )
      {
        independent_to_result[ itr->GetIndependentDatasetRetriever()].PushBack( itr->GetResult());
      }

      // determine the directionality of improvement to create
      opti::ImprovementTypeEnum
        obj_improvement_type( cv_info.FirstElement().GetImprovementType());

      const size_t results_size( DATASET.GetResultSize());
      storage::Vector< linal::Matrix< char> > prediction_classifications
      (
        n_models,
        linal::Matrix< char>( size_t( 1), results_size, '\0')
      );

      // determine the cutoff for each independent set
      for
      (
        storage::Map< util::ObjectDataLabel, storage::Vector< float> >::iterator
          itr( independent_to_result.Begin()), itr_end( independent_to_result.End());
        itr != itr_end;
        ++itr
      )
      {
        itr->second.Sort( std::less< float>());
        const size_t size( itr->second.GetSize());
        const size_t median( size / 2);
        const float low_value( itr->second.FirstElement());
        const float median_value( itr->second( median));
        const float high_value( itr->second.LastElement());
        if( itr->second.FirstElement() == itr->second.LastElement())
        {
          independent_to_median[ itr->first] = util::GetUndefined< float>();
        }
        else if( median_value == low_value || high_value - median_value > median_value - low_value)
        {
          size_t new_median( median + 1);
          while( median_value == itr->second( new_median))
          {
            ++new_median;
          }
          independent_to_median[ itr->first] = 0.5 * ( median_value + itr->second( new_median));
        }
        else
        {
          size_t new_median( median - 1);
          while( median_value == itr->second( new_median))
          {
            --new_median;
          }
          independent_to_median[ itr->first] = 0.5 * ( median_value + itr->second( new_median));
        }
      }

      // create a vector with size-t's of 0 for poorly-performing models in the CV, 1 for better-performing models
      // undefined for models for which this information is not available
      storage::Vector< std::string> model_score( n_models, std::string( results_size, '\0'));
      storage::List< CrossValidationInfo>::const_iterator itr_cv( cv_info.Begin());
      bool larger_is_better( obj_improvement_type == opti::e_LargerEqualIsBetter);
      for( size_t model_number( 0); model_number < n_models; ++model_number, ++itr_cv)
      {
        // get the independent result median
        const float ind_result_median( independent_to_median[ itr_cv->GetIndependentDatasetRetriever()]);
        if( !util::IsDefined( ind_result_median))
        {
          // undefined median, continue
          continue;
        }
        if( larger_is_better == ( itr_cv->GetResult() > ind_result_median))
        {
          prediction_classifications( model_number) = 'P';
          std::fill( model_score( model_number).begin(), model_score( model_number).end(), 'P');
        }
        else
        {
          prediction_classifications( model_number) = 'n';
          std::fill( model_score( model_number).begin(), model_score( model_number).end(), 'n');
        }
      }
      m_Scorer.InitializeBalancing( prediction_classifications, feature_size);
      return m_Scorer.Score( weight_effects, model_score);
    }

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERR_STREAM stream to write errors out to
    bool ScoreDatasetNeuralNetworkWeights::ReadInitializerSuccessHook
    (
      const util::ObjectDataLabel &LABEL,
      std::ostream &ERR_STREAM
    )
    {
      // get the models if they are not already in memory
      storage::Vector< util::OwnPtr< NeuralNetwork> > &models( s_Models[ GetString()]);

      if( m_Key.empty())
      {
        // get all the keys for this storage
        storage::Vector< std::string> keys( m_ModelStorage->GetAllKeys());

        BCL_MessageVrb
        (
          "Found " + util::Format()( keys.GetSize()) + " models for neural network weights"
        );

        if( models.IsEmpty() && !keys.IsEmpty())
        {
          // no models already loaded
          RetrieveInterface::t_Container interfaces( m_ModelStorage->RetrieveEnsemble( keys));
          models = storage::Vector< util::OwnPtr< NeuralNetwork> >( interfaces.Begin(), interfaces.End());
        }
        if( keys.GetSize() != models.GetSize())
        {
          ERR_STREAM << "# of neural networks in " << m_ModelStorage->GetString() << " changed; aborting";
          return false;
        }
        else if( keys.IsEmpty())
        {
          ERR_STREAM << "No neural networks found in storage! " << m_ModelStorage->GetString();
          return false;
        }
      }
      else if( models.IsEmpty())
      {
        // no models already loaded
        models.PushBack( m_ModelStorage->Retrieve( m_Key));
        if( !models.LastElement().IsDefined())
        {
          ERR_STREAM << m_Key << " is not a key for " << m_ModelStorage->GetString();
          return false;
        }
      }

      m_Models = util::ToSiPtr( models);
      return true;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer ScoreDatasetNeuralNetworkWeights::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription( "Uses neural network weights as a derivative proxy for scoring");

      parameters.AddInitializer
      (
        "storage",
        "type of storage for models",
        io::Serialization::GetAgent( &m_ModelStorage)
      );

      parameters.AddInitializer
      (
        "key",
        "optional key of the desired model in the storage",
        io::Serialization::GetAgent( &m_Key),
        ""
      );

      parameters.AddInitializer
      (
        "weights",
        "Change the weighting of various measures of the neural networks' weights",
        io::Serialization::GetAgent( &m_Scorer)
      );
      return parameters;
    }

  } // namespace model
} // namespace bcl

