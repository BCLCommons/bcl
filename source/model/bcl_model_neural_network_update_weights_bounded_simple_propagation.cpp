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
#include "model/bcl_model_neural_network_update_weights_bounded_simple_propagation.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {

    // add the interface to the set of known implementations
    const util::SiPtr< const util::ObjectInterface> NeuralNetworkUpdateWeightsBoundedSimplePropagation::s_Instance
    (
      util::Enumerated< NeuralNetworkUpdateWeightsInterface>::AddInstance
      (
        new NeuralNetworkUpdateWeightsBoundedSimplePropagation()
      )
    );

    //! @brief default constructor
    NeuralNetworkUpdateWeightsBoundedSimplePropagation::NeuralNetworkUpdateWeightsBoundedSimplePropagation() :
      m_Alpha( 0.5),
      m_Eta( 10.0),
      m_Min( math::GetLowestBoundedValue< float>()),
      m_Max( math::GetHighestBoundedValue< float>())
    {
    }

    //! @brief copy constructor
    //! @return a new NeuralNetworkUpdateWeightsBoundedSimplePropagation copied from this instance
    NeuralNetworkUpdateWeightsBoundedSimplePropagation *NeuralNetworkUpdateWeightsBoundedSimplePropagation::Clone() const
    {
      return new NeuralNetworkUpdateWeightsBoundedSimplePropagation( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &NeuralNetworkUpdateWeightsBoundedSimplePropagation::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &NeuralNetworkUpdateWeightsBoundedSimplePropagation::GetAlias() const
    {
      static const std::string s_Name( "BoundedSimple");
      return s_Name;
    }

    //! @brief initialize the update weights; should initialize any internal data structures to the right size
    //! @param SIZE the number of weights/changes/slopes/prevslopes to update
    void NeuralNetworkUpdateWeightsBoundedSimplePropagation::Initialize( const size_t &SIZE)
    {
      m_ChangeSlopes = linal::Vector< float>( SIZE, float( 0.0));
    }

    //! @brief Set the changes array; this is intended solely for testing purposes
    //! @param CHANGES the new changes array
    void NeuralNetworkUpdateWeightsBoundedSimplePropagation::SetChanges( const linal::VectorConstInterface< float> &CHANGES)
    {
      m_ChangeSlopes = CHANGES;
    }

    //! @brief Get the changes array; this is intended solely for testing purposes
    //! @return the changes array
    const linal::Vector< float> &NeuralNetworkUpdateWeightsBoundedSimplePropagation::GetChanges() const
    {
      return m_ChangeSlopes;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief this is a function used internally to update the weights of a neural network
    //! @param WEIGHTS a vector iterator (const pointer) to the first weight to be updated
    //! @param SLOPES  a vector iterator (const pointer) to the first slope to be updated
    void NeuralNetworkUpdateWeightsBoundedSimplePropagation::operator ()
    (
      float *const &WEIGHTS,
      float *const &SLOPES
    )
    {
      const size_t number_connection_weights( m_ChangeSlopes.GetSize());
      float *const changes( m_ChangeSlopes.Begin());
      // simple propagation algorithm
      for( size_t i( 0); i < number_connection_weights; ++i)
      {
        if( util::IsNaN( SLOPES[ i]))
        {
          continue;
        }

        changes[ i] = m_Eta * SLOPES[ i] + m_Alpha * changes[ i];
        WEIGHTS[ i] += changes[ i];
        SLOPES[ i] = 0;
        if( WEIGHTS[ i] > m_Max)
        {
          WEIGHTS[ i] = m_Max;
        }
        else if( WEIGHTS[ i] < m_Min)
        {
          WEIGHTS[ i] = m_Min;
        }
      }
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer NeuralNetworkUpdateWeightsBoundedSimplePropagation::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "uses simple backpropagation (see http://en.wikipedia.org/wiki/Backpropagation)"
      );

      parameters.AddInitializer
      (
        "eta",
        "learning rate",
        io::Serialization::GetAgentWithRange( &m_Eta, 0.0, 1.0),
        "0.1"
      );
      parameters.AddInitializer
      (
        "alpha",
        "momentum; how long a change in weights persists",
        io::Serialization::GetAgentWithRange( &m_Alpha, 0.0, 1.0),
        "0.5"
      );
      parameters.AddInitializer
      (
        "min",
        "minimum value for the weight; if it goes below this, it will be reset to this",
        io::Serialization::GetAgent( &m_Min),
        util::Format()( math::GetLowestBoundedValue< float>())
      );
      parameters.AddInitializer
      (
        "max",
        "maximum value for the weight; if it goes above this, it will be reset to this",
        io::Serialization::GetAgent( &m_Max),
        util::Format()( math::GetHighestBoundedValue< float>())
      );
      return parameters;
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &NeuralNetworkUpdateWeightsBoundedSimplePropagation::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Alpha, ISTREAM);
      io::Serialize::Read( m_Eta, ISTREAM);
      io::Serialize::Read( m_ChangeSlopes, ISTREAM);
      io::Serialize::Read( m_Min, ISTREAM);
      io::Serialize::Read( m_Max, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &NeuralNetworkUpdateWeightsBoundedSimplePropagation::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Alpha, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Eta, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_ChangeSlopes, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Min, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Max, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  } // namespace model
} // namespace bcl
