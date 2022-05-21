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
#include "model/bcl_model_neural_network_update_weights_resilient_propagation.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {

    // add the interface to the set of known implementations
    const util::SiPtr< const util::ObjectInterface> NeuralNetworkUpdateWeightsResilientPropagation::s_Instance
    (
      util::Enumerated< NeuralNetworkUpdateWeightsInterface>::AddInstance
      (
        new NeuralNetworkUpdateWeightsResilientPropagation()
      )
    );

    //! @brief default constructor
    NeuralNetworkUpdateWeightsResilientPropagation::NeuralNetworkUpdateWeightsResilientPropagation() :
      m_PreviousSlopes(),
      m_MaxWeightChange( 50.0),
      m_MinWeightChange( 0.001)
    {
    }

    //! @brief copy constructor
    //! @return a new NeuralNetworkUpdateWeightsResilientPropagation copied from this instance
    NeuralNetworkUpdateWeightsResilientPropagation *NeuralNetworkUpdateWeightsResilientPropagation::Clone() const
    {
      return new NeuralNetworkUpdateWeightsResilientPropagation( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &NeuralNetworkUpdateWeightsResilientPropagation::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &NeuralNetworkUpdateWeightsResilientPropagation::GetAlias() const
    {
      static const std::string s_Name( "Resilient");
      return s_Name;
    }

    //! @brief initialize the update weights; should initialize any internal data structures to the right size
    //! @param SIZE the number of weights/changes/slopes/prevslopes to update
    void NeuralNetworkUpdateWeightsResilientPropagation::Initialize( const size_t &SIZE)
    {
      // set the previous slopes to 0.0125
      m_PreviousSlopes = linal::Vector< float>( SIZE, 0.0125);
      m_ChangeSlopes = linal::Vector< float>( SIZE, 0.0);
    } // Initialize

    //! @brief Set the changes array; this is intended solely for testing purposes
    //! @param CHANGES the new changes array
    void NeuralNetworkUpdateWeightsResilientPropagation::SetChanges( const linal::VectorConstInterface< float> &CHANGES)
    {
      m_ChangeSlopes = CHANGES;
    }

    //! @brief Get the changes array; this is intended solely for testing purposes
    //! @return the changes array
    const linal::Vector< float> &NeuralNetworkUpdateWeightsResilientPropagation::GetChanges() const
    {
      return m_ChangeSlopes;
    }

    //! @brief this is a function used internally to update the weights of a neural network
    //! @param WEIGHTS a vector iterator (const pointer) to the first weight to be updated
    //! @param SLOPES  a vector iterator (const pointer) to the first slope to be updated
    void NeuralNetworkUpdateWeightsResilientPropagation::operator ()( float *const &WEIGHTS, float *const &SLOPES)
    {
      const float increase_factor( 1.2);
      const float decrease_factor( 0.5);

      float *const changes( m_ChangeSlopes.Begin());
      float *const previous_slopes( m_PreviousSlopes.Begin());

      // computes new change
      for( size_t index( 0), size( m_PreviousSlopes.GetSize()); index < size; ++index)
      {
        float &slope( SLOPES[ index]);
        if( util::IsNaN( slope))
        {
          continue;
        }

        float &weight( WEIGHTS[ index]);
        float &change( changes[ index]);
        float &previous_slope( previous_slopes[ index]);
        const float same_sign( previous_slope * slope);

        // this slope and previous slope have the same sign, increase the change accordingly
        if( same_sign > std::numeric_limits< float>::min())
        {
          change = std::min( std::max( change, m_MinWeightChange) * increase_factor, m_MaxWeightChange);
          weight += slope < float( 0.0) ? -change : change;
        }
        // this slope and previous slope have different signs
        else if( same_sign < -std::numeric_limits< float>::min())
        {
          // change may not be zero because then the training will stop
          change = std::max( change, m_MinWeightChange) * decrease_factor;
          slope = 0.0;
        }
        // either the current slope or the previous slope was zero
        // usually its the previous slope that is zero, which happens if, in the last iteration, slope and
        // previous slope had opposite signs
        else
        {
          weight += slope < float( 0.0) ? -change : change;
        }

        // the previous slope for the next round
        previous_slope = slope;
      }
    } // ResilientUpdateWeights

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer NeuralNetworkUpdateWeightsResilientPropagation::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "uses the resilient propagation (rprop) algorithm (see http://en.wikipedia.org/wiki/Rprop)"
      );
      parameters.AddInitializer
      (
        "max change",
        "max weight change allowed per iteration",
        io::Serialization::GetAgent( &m_MaxWeightChange),
        "50.0"
      );
      parameters.AddInitializer
      (
        "min change",
        "min weight change allowed per iteration",
        io::Serialization::GetAgent( &m_MinWeightChange),
        "0.001"
      );
      return parameters;
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &NeuralNetworkUpdateWeightsResilientPropagation::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_PreviousSlopes, ISTREAM);
      io::Serialize::Read( m_ChangeSlopes, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &NeuralNetworkUpdateWeightsResilientPropagation::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_PreviousSlopes, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_ChangeSlopes, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  } // namespace model
} // namespace bcl
