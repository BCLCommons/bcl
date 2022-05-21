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

#ifndef BCL_MODEL_NEURAL_NETWORK_PERTURBATION_INTERFACE_H_
#define BCL_MODEL_NEURAL_NETWORK_PERTURBATION_INTERFACE_H_

// include the namespace header
#include "bcl_model.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_matrix.h"
#include "util/bcl_util_enumerated.h"
#include "util/bcl_util_serializable_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class NeuralNetworkPerturbationInterface
    //! @brief an abstracted method of perturbing a neural network, specifically its weights
    //!
    //! @remarks example unnecessary
    //! @author mendenjl
    //! @date May 03, 2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API NeuralNetworkPerturbationInterface :
      public util::SerializableInterface
    {

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief virtual copy constructor
      //! @return a new NeuralNetworkPerturbationInterface copied from this instance
      virtual NeuralNetworkPerturbationInterface *Clone() const = 0;

    ///////////////
    // operators //
    ///////////////

      //! @brief this is a function used internally to update the weights of a neural network
      //! @param WEIGHTS a vector iterator (const pointer) to the first weight to be updated
      virtual void operator ()( linal::Matrix< float> &WEIGHTS) = 0;

    }; // class NeuralNetworkPerturbationInterface

  } // namespace model
} // namespace bcl

#endif // BCL_MODEL_NEURAL_NETWORK_PERTURBATION_INTERFACE_H_
