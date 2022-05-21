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
#include "model/bcl_model_neural_network_perturb_attenuate.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {

    // add the interface to the set of known implementations
    const util::SiPtr< const util::ObjectInterface> NeuralNetworkPerturbAttenuate::s_Instance
    (
      util::Enumerated< NeuralNetworkPerturbationInterface>::AddInstance
      (
        new NeuralNetworkPerturbAttenuate()
      )
    );

    //! @brief copy constructor
    //! @return a new NeuralNetworkPerturbAttenuate copied from this instance
    NeuralNetworkPerturbAttenuate *NeuralNetworkPerturbAttenuate::Clone() const
    {
      return new NeuralNetworkPerturbAttenuate( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &NeuralNetworkPerturbAttenuate::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &NeuralNetworkPerturbAttenuate::GetAlias() const
    {
      static const std::string s_Name( "Attenuate");
      return s_Name;
    }

    //! @brief this is a function used internally to update the weights of a neural network
    //! @param WEIGHTS a matrix of weights
    void NeuralNetworkPerturbAttenuate::operator ()( linal::Matrix< float> &WEIGHTS)
    {
      // attenuate the weights
      for( linal::Matrix< float>::iterator itr( WEIGHTS.Begin()), itr_end( WEIGHTS.End()); itr != itr_end; ++itr)
      {
        // get the weight
        float &weight( *itr);

        // compute absolute weight
        const float abs_weight( math::Absolute( weight));

        // determine by how much to attenuate the weight
        const float weight_attenuation( m_Attenuation( abs_weight));

        // if the attenuation would be larger than the current weight, set the weight to zero
        if( weight_attenuation > abs_weight)
        {
          weight = 0.0;
        }
        // otherwise, add or subtract the weight attenuation to move the weight closer to 0
        else if( weight < float( 0.0))
        {
          weight += weight_attenuation;
        }
        else
        {
          weight -= weight_attenuation;
        }
      }
    } // operator ()

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer NeuralNetworkPerturbAttenuate::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Attenuates weights directly according to the absolute value of the weight"
      );
      parameters.AddInitializer
      (
        "",
        "Coefficients of polynomial used to attenuate weights, in ascending order of degree",
        io::Serialization::GetAgent( &m_Attenuation)
      );
      return parameters;
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &NeuralNetworkPerturbAttenuate::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Attenuation, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &NeuralNetworkPerturbAttenuate::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Attenuation, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  } // namespace model
} // namespace bcl
