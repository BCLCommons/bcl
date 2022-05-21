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

#ifndef BCL_MODEL_PRETRAIN_NEURAL_NETWORK_FROM_FILE_H_
#define BCL_MODEL_PRETRAIN_NEURAL_NETWORK_FROM_FILE_H_

// include the namespace header
#include "bcl_model.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_model_pretrain_neural_network_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class PretrainNeuralNetworkFromFile
    //! @brief Iterate that was read in from a file
    //!
    //! @see @link example_model_pretrain_neural_network_from_file.cpp @endlink
    //! @author mendenjl
    //! @date Aug 14, 2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API PretrainNeuralNetworkFromFile :
      public PretrainNeuralNetworkInterface
    {

    private:

    //////////
    // data //
    //////////

      // filename for neural network
      std::string m_FileName;

      // pointer to neural network implementation
      util::ShPtr< NeuralNetwork> m_NeuralNetwork;

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief Clone function
      //! @return pointer to new PretrainNeuralNetworkFromFile
      PretrainNeuralNetworkFromFile *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief get a short name for this class
      //! @return a short name for this class
      const std::string &GetAlias() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief the main operation, pretrains a neural network
      //! @param DATA the data for use in pretraining
      //! @param OBJECTIVE ShPtr to the objective function for the network
      util::ShPtr< NeuralNetwork> PretrainNetwork
      (
        util::ShPtr< descriptor::Dataset> &DATA,
        const util::ShPtr< ObjectiveFunctionWrapper> &OBJECTIVE
      );

    ///////////////
    // operators //
    ///////////////

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      //! @brief Set the members of this property from the given LABEL
      //! @param LABEL the label to parse
      //! @param ERROR_STREAM the stream to write errors to
      bool ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM);

    }; // class PretrainNeuralNetworkFromFile

  } // namespace model
} // namespace bcl

#endif // BCL_MODEL_PRETRAIN_NEURAL_NETWORK_FROM_FILE_H_
