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
#include "model/bcl_model_neural_network_perturb_max_norm.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"
#include "linal/bcl_linal_vector_operations.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {

    // add the interface to the set of known implementations
    const util::SiPtr< const util::ObjectInterface> NeuralNetworkPerturbMaxNorm::s_Instance
    (
      util::Enumerated< NeuralNetworkPerturbationInterface>::AddInstance
      (
        new NeuralNetworkPerturbMaxNorm()
      )
    );

    //! @brief copy constructor
    //! @return a new NeuralNetworkPerturbMaxNorm copied from this instance
    NeuralNetworkPerturbMaxNorm *NeuralNetworkPerturbMaxNorm::Clone() const
    {
      return new NeuralNetworkPerturbMaxNorm( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &NeuralNetworkPerturbMaxNorm::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &NeuralNetworkPerturbMaxNorm::GetAlias() const
    {
      static const std::string s_Name( "MaxNorm");
      return s_Name;
    }

    //! @brief this is a function used internally to update the weights of a neural network
    //! @param WEIGHTS a matrix of weights
    void NeuralNetworkPerturbMaxNorm::operator ()( linal::Matrix< float> &WEIGHTS)
    {
      // perform normalization on inputs 1st: normalize all rows such that L2(row) <= m_MaxNormIncoming

      // perform the max-norm operation on the rows, simultaneously computing the column-wise sum of squares
      linal::Vector< float> col_norms( WEIGHTS.GetNumberCols(), float( 0.0));
      size_t n_inputs_normalized( 0);
      for( size_t row( 0), n_rows( WEIGHTS.GetNumberRows()); row < n_rows; ++row)
      {
        // get a reference to the row
        linal::VectorReference< float> row_ref( WEIGHTS.GetRow( row));

        // compute norm
        const float norm( row_ref.Norm());

        // if norm > m_MaxNormIncoming, scale
        if( norm > m_MaxNormIncoming)
        {
          row_ref *= m_MaxNormIncoming / norm;
          ++n_inputs_normalized;
        }

        col_norms += linal::SqrVector( row_ref);
      }

      // compute column-wise norms by taking the sqrt of col_norms
      // Then, transform col_norms into the necessary multipling vector for each row in the matrix
      size_t n_outputs_normalized( 0);
      for( size_t col( 0), n_cols( col_norms.GetSize()); col < n_cols; ++col)
      {
        const float col_norm( math::Sqrt( col_norms( col)));
        if( col_norm > m_MaxNormOutgoing)
        {
          col_norms( col) = m_MaxNormOutgoing / col_norm;
          ++n_outputs_normalized;
        }
        else
        {
          // norm is sufficiently small, set to 1 (no multiplication)
          col_norms( col) = 1.0;
        }
      }

      if( n_inputs_normalized + n_outputs_normalized)
      {
        BCL_MessageStd
        (
          "# of inputs normalized: " + util::Format()( n_inputs_normalized)
          + "# of outputs normalized: " + util::Format()( n_outputs_normalized)
        );
      }

      // check if there are any columns to normalize
      if( !n_outputs_normalized)
      {
        return;
      }

      // normalize by multiplying each row of WEIGHTS by col_norms
      for( size_t row( 0), n_rows( WEIGHTS.GetNumberRows()); row < n_rows; ++row)
      {
        // get a reference to the row
        linal::VectorReference< float> row_ref( WEIGHTS.GetRow( row));

        // multiply each row by col_norms
        linal::ElementwiseMultiply( row_ref, col_norms);
      }
    } // operator ()

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer NeuralNetworkPerturbMaxNorm::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Limits the incoming weight norm (L2) for the neuron; if weights exceed this value, all values are scaled"
      );
      parameters.AddInitializer
      (
        "in",
        "Max euclidean norm for hidden neuron input weights",
        io::Serialization::GetAgent( &m_MaxNormIncoming)
      );
      parameters.AddInitializer
      (
        "out",
        "Max euclidean norm for input neuron output to the hidden layer (normally substantially larger than in)",
        io::Serialization::GetAgent( &m_MaxNormOutgoing),
        "10"
      );
      return parameters;
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &NeuralNetworkPerturbMaxNorm::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_MaxNormIncoming, ISTREAM);
      io::Serialize::Read( m_MaxNormOutgoing, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &NeuralNetworkPerturbMaxNorm::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_MaxNormIncoming, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_MaxNormOutgoing, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  } // namespace model
} // namespace bcl
