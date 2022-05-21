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

#ifndef BCL_MODEL_NEURAL_NETWORK_UPDATE_WEIGHTS_INTERFACE_H_
#define BCL_MODEL_NEURAL_NETWORK_UPDATE_WEIGHTS_INTERFACE_H_

// include the namespace header
#include "bcl_model.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_vector.h"
#include "util/bcl_util_enumerated.h"
#include "util/bcl_util_serializable_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class NeuralNetworkUpdateWeightsInterface
    //! @brief an abstracted method of updating a neural networks weights
    //!
    //! @remarks example unnecessary
    //! @author mendenjl
    //! @date May 12, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API NeuralNetworkUpdateWeightsInterface :
      public util::SerializableInterface
    {

    public:

    //////////
    // data //
    //////////

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief virtual copy constructor
      //! @return a new NeuralNetworkUpdateWeightsInterface copied from this instance
      virtual NeuralNetworkUpdateWeightsInterface *Clone() const = 0;

    /////////////////
    // data access //
    /////////////////

    ////////////////
    // operations //
    ////////////////

      //! @brief initialize the update weights; should initialize any internal data structures to the right size
      //! @param SIZE the number of weights/changes/slopes/prevslopes to update
      virtual void Initialize( const size_t &SIZE) = 0;

      //! @brief Get the changes array; this is intended solely for testing purposes
      //! @return the changes array
      virtual const linal::Vector< float> &GetChanges() const = 0;

      //! @brief Set the changes array; this is intended solely for testing purposes
      //! @param CHANGES the new changes array
      virtual void SetChanges( const linal::VectorConstInterface< float> &CHANGES) = 0;

    ///////////////
    // operators //
    ///////////////

      //! @brief this is a function used internally to update the weights of a neural network
      //! @param WEIGHTS a vector iterator (const pointer) to the first weight to be updated
      //! @param SLOPES  a vector iterator (const pointer) to the first slope to be updated
      virtual void operator ()( float *const &WEIGHTS, float *const &SLOPES) = 0;

    }; // class NeuralNetworkUpdateWeightsInterface

  } // namespace model
} // namespace bcl

#endif // BCL_MODEL_NEURAL_NETWORK_UPDATE_WEIGHTS_INTERFACE_H_
