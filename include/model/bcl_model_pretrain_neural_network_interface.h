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

#ifndef BCL_MODEL_PRETRAIN_NEURAL_NETWORK_INTERFACE_H_
#define BCL_MODEL_PRETRAIN_NEURAL_NETWORK_INTERFACE_H_

// include the namespace header
#include "bcl_model.h"

// include other forward headers - sorted alphabetically
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_model_has_labels_base.h"
#include "bcl_model_neural_network.h"
#include "descriptor/bcl_descriptor_dataset.h"
#include "util/bcl_util_enumerated.h"
#include "util/bcl_util_serializable_interface.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class PretrainNeuralNetworkInterface
    //! @brief an interface for objects that have some way of pretraining an initial neural network
    //!
    //! @see @link example_model_pretrain_neural_network_interface.cpp @endlink
    //! @author mendenjl
    //! @date Aug 14, 2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API PretrainNeuralNetworkInterface :
      public util::SerializableInterface,
      virtual public HasLabelsBase
    {

    public:

    /////////////////
    // data access //
    /////////////////

      //! @brief virtual copy constructor
      virtual PretrainNeuralNetworkInterface *Clone() const = 0;

    ////////////////
    // operations //
    ////////////////

      //! @brief the main operation, pretrains a neural network
      //! @param DATA the data for use in pretraining
      //! @param OBJECTIVE ShPtr to the objective function for the network
      virtual util::ShPtr< NeuralNetwork> PretrainNetwork
      (
        util::ShPtr< descriptor::Dataset> &DATA,
        const util::ShPtr< ObjectiveFunctionWrapper> &OBJECTIVE
      ) = 0;

    }; // class PretrainNeuralNetworkInterface

  } // namespace model
} // namespace bcl

#endif // BCL_MODEL_PRETRAIN_NEURAL_NETWORK_INTERFACE_H_
