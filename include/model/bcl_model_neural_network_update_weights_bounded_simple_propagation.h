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

#ifndef BCL_MODEL_NEURAL_NETWORK_UPDATE_WEIGHTS_BOUNDED_SIMPLE_PROPAGATION_H_
#define BCL_MODEL_NEURAL_NETWORK_UPDATE_WEIGHTS_BOUNDED_SIMPLE_PROPAGATION_H_

// include the namespace header
#include "bcl_model.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_model_neural_network_update_weights_interface.h"
#include "linal/bcl_linal_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class NeuralNetworkUpdateWeightsBoundedSimplePropagation
    //! @brief provides the weight update function for training a neural network with the
    //!        simple propagation update algorithm
    //!
    //! @see @link example_model_neural_network_update_weights_bounded_simple_propagation.cpp @endlink
    //! @author mendenjl
    //! @date Jun 27, 2017
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API NeuralNetworkUpdateWeightsBoundedSimplePropagation :
      public NeuralNetworkUpdateWeightsInterface
    {
    private:

    //////////
    // data //
    //////////

      //! learning momentum
      float m_Alpha;

      //! learning rate
      float m_Eta;

      //! Minimum value; if the value of weight drops below this, it is automatically reset to this
      float m_Min;

      //! Maximum value; if the weight climbs above this, it is automatically reset to this
      float m_Max;

      linal::Vector< float> m_ChangeSlopes;   //!< Changes to previous slopes

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      NeuralNetworkUpdateWeightsBoundedSimplePropagation();

      //! @brief copy constructor
      //! @return a new NeuralNetworkUpdateWeightsBoundedSimplePropagation copied from this instance
      NeuralNetworkUpdateWeightsBoundedSimplePropagation *Clone() const;

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

      //! @brief initialize the update weights; should initialize any internal data structures to the right size
      //! @param SIZE the number of weights/changes/slopes/prevslopes to update
      void Initialize( const size_t &SIZE);

      //! @brief Get the changes array; this is intended solely for testing purposes
      //! @return the changes array
      const linal::Vector< float> &GetChanges() const;

      //! @brief Set the changes array; this is intended solely for testing purposes
      //! @param CHANGES the new changes array
      void SetChanges( const linal::VectorConstInterface< float> &CHANGES);

    ///////////////
    // operators //
    ///////////////

      //! @brief this is a function used internally to update the weights of a neural network
      //! @param WEIGHTS a vector iterator (const pointer) to the first weight to be updated
      //! @param SLOPES  a vector iterator (const pointer) to the first slope to be updated
      void operator()( float *const &WEIGHTS, float *const &SLOPES);

    //////////////////////
    // input and output //
    //////////////////////

    protected:

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class IterateSimpleProagation

  } // namespace model
} // namespace bcl

#endif // BCL_MODEL_NEURAL_NETWORK_UPDATE_WEIGHTS_BOUNDED_SIMPLE_PROPAGATION_H_
