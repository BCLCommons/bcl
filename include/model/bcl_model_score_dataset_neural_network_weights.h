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

#ifndef BCL_MODEL_SCORE_DATASET_NEURAL_NETWORK_WEIGHTS_H
#define BCL_MODEL_SCORE_DATASET_NEURAL_NETWORK_WEIGHTS_H

// include the namespace header
#include "bcl_model.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_model_retrieve_interface.h"
#include "bcl_model_score_dataset_interface.h"
#include "bcl_model_score_derivative_ensemble.h"
#include "linal/bcl_linal_vector.h"
#include "storage/bcl_storage_map.h"
#include "util/bcl_util_implementation.h"
#include "util/bcl_util_sh_ptr_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ScoreDatasetNeuralNetworkWeights
    //! @brief calculates the saliency of all input features by taking the sum of absolute weights emanating from each input neuron
    //!
    //! @author mendenjl
    //! @see @link example_model_score_dataset_neural_network_weights.cpp @endlink
    //! @date May 03, 2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ScoreDatasetNeuralNetworkWeights :
      public ScoreDatasetInterface
    {
    private:

    //////////
    // data //
    //////////

      //! Model storage; to retrieve models and descriptors
      util::Implementation< RetrieveInterface> m_ModelStorage;

      //! Optional key of the desired model
      std::string m_Key;

      //! scorer for the nominal derivatives calculated by this function
      ScoreDerivativeEnsemble m_Scorer;

      //! pointer to vector of models used by this property
      util::SiPtr< const util::ShPtrVector< NeuralNetwork> > m_Models;

      //! Map from initializer string/type to the shptrvector of model interfaces; saves loading potentially
      //! gigantic models repetitively
      static storage::Map< std::string, storage::Vector< util::OwnPtr< NeuralNetwork> > > s_Models;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief Clone function
      //! @return pointer to new ScoreDatasetNeuralNetworkWeights
      ScoreDatasetNeuralNetworkWeights *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief get a short name for this class
      //! @return a short name for this class
      const std::string &GetAlias() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief score a given dataset
      //! @param DATASET dataset of interest
      //! @return scores of the dataset
      linal::Vector< float> Score( const descriptor::Dataset &DATASET) const;

    ///////////////
    // operators //
    ///////////////

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief Set the members of this property from the given LABEL
      //! @param LABEL the label to parse
      //! @param ERR_STREAM stream to write errors out to
      bool ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERR_STREAM);

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    }; // class ScoreDatasetNeuralNetworkWeights

  } // namespace model
} // namespace bcl

#endif // BCL_MODEL_SCORE_DATASET_NEURAL_NETWORK_WEIGHTS_H

