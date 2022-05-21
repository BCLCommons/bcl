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
#include "model/bcl_model_support_vector_machine_multi_output.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"
#include "linal/bcl_linal_vector_operations.h"
#include "model/bcl_model_rescale_feature_data_set.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {

  //////////
  // data //
  //////////

    //! the default input range for SVM normalization range
    const math::Range< float> SupportVectorMachineMultiOutput::s_DefaultInputRange( 0, 1);

    // add the interface to the set of known implementations
    const util::SiPtr< const util::ObjectInterface> SupportVectorMachineMultiOutput::s_Instance
    (
      GetObjectInstances().AddInstance( new SupportVectorMachineMultiOutput())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief standard constructor
    SupportVectorMachineMultiOutput::SupportVectorMachineMultiOutput() :
      m_Beta(),
      m_Bias( linal::Vector< float>( 1, float(0))),
      m_Kernel(),
      m_NumberBoundSupportVectors(),
      m_NumberSupportVectors(),
      m_RescaleInput(),
      m_RescaleOutput(),
      m_SupportVectors()
    {
    }

    //! @brief standard constructor
    //! @param KERNEL used Kernel Function
    //! @param RESCALE_INPUT Rescale Function for inputs (descriptors)
    //! @param RESCALE_OUTPUT Rescale Function for outputs (results)
    SupportVectorMachineMultiOutput::SupportVectorMachineMultiOutput
    (
      const util::Implementation< SupportVectorKernelBase> &KERNEL,
      const util::ShPtr< RescaleFeatureDataSet> &RESCALE_INPUT,
      const util::ShPtr< RescaleFeatureDataSet> &RESCALE_OUTPUT
    ) :
      m_Beta( storage::Vector< linal::Vector< float> >()),
      m_Bias( linal::Vector< float>( 1, float( 0))),
      m_Kernel( KERNEL),
      m_NumberBoundSupportVectors( size_t( 0)),
      m_NumberSupportVectors( size_t( 0)),
      m_RescaleInput( RESCALE_INPUT),
      m_RescaleOutput( RESCALE_OUTPUT),
      m_SupportVectors( FeatureDataSet< float>( linal::MatrixReference< float>()))
    {
    }

    //! @brief standard constructor
    //! @param KERNEL used Kernel Function
    SupportVectorMachineMultiOutput::SupportVectorMachineMultiOutput
    (
      const util::Implementation< SupportVectorKernelBase> &KERNEL
    ) :
      m_Beta(),
      m_Bias( linal::Vector< float>( 1, float(0))),
      m_Kernel( KERNEL),
      m_NumberBoundSupportVectors( size_t( 0)),
      m_NumberSupportVectors( size_t( 0)),
      m_RescaleInput( new RescaleFeatureDataSet()),
      m_RescaleOutput( new RescaleFeatureDataSet()),
      m_SupportVectors( FeatureDataSet< float>())
    {
    }

    //! @brief virtual copy constructor
    SupportVectorMachineMultiOutput *SupportVectorMachineMultiOutput::Clone() const
    {
      return new SupportVectorMachineMultiOutput( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &SupportVectorMachineMultiOutput::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &SupportVectorMachineMultiOutput::GetAlias() const
    {
      static const std::string s_Name( "SupportVectorMachineMultiOutput");
      return s_Name;
    }

    //! @brief Get a description for what this class does (used when writing help)
    //! @return a description for what this class does (used when writing help)
    const std::string &SupportVectorMachineMultiOutput::GetClassDescription() const
    {
      static const std::string s_Description( "see http://en.wikipedia.org/wiki/Support_vector_machine");
      return s_Description;
    }

  //////////////
  // operator //
  //////////////

  ////////////////
  // operations //
  ////////////////

    //! @brief Set the scaling of a feature set according to the model
    //! @param FEATURES feature set of interest
    //! @note this allows for external classes that own a dataset to ensure that a new dataset is never created
    //!       when operator() is called
    void SupportVectorMachineMultiOutput::Rescale( FeatureDataSet< float> &FEATURE) const
    {
      if( !FEATURE.IsRescaled() || *FEATURE.GetScaling() != *m_RescaleInput)
      {
        FEATURE.DeScale();
        FEATURE.Rescale( *m_RescaleInput);
      }
    }

    //! @brief predict result with model using a rescaled feature vector
    //! @param FEATURE normalized or rescaled feature vector
    //! @return predcited normalized result vector using a model
    FeatureDataSet< float> SupportVectorMachineMultiOutput::PredictWithoutRescaling
    (
      const FeatureDataSetInterface< float> &FEATURE
    ) const
    {
      const size_t num_points( FEATURE.GetNumberFeatures());
      linal::Matrix< float> prediction( num_points, m_Bias.GetSize(), float( 0.0));

      for( size_t i( 0); i < num_points; ++i)
      {
        // initialize classification result with bias
        linal::Vector< float> classification( m_Bias);

        // progress counter for accessing data set with support vectors
        size_t progress( 0);

        //iterator over all LaGrange multiplier alpha
        for
        (
          storage::Vector< linal::Vector< float> >::const_iterator iter_beta( GetBeta().Begin()),
          iter_beta_end( GetBeta().End());
          iter_beta != iter_beta_end;
          ++iter_beta, ++progress
        )
        {
          // calculate classification value of beta and kernel value of unknown vector and a support vector
          classification +=
            *iter_beta * m_Kernel->operator()( FEATURE( i), m_SupportVectors( progress));
        }
        prediction.ReplaceRow( i, classification);
      }

      // classification result
      return FeatureDataSet< float>( prediction);

    }

    //! @brief predict result with model using a NOT rescaled feature vector
    //! @param FEATURE not rescaled feature vector
    //! @return predcited result vector using a model
    FeatureDataSet< float> SupportVectorMachineMultiOutput::operator()( const FeatureDataSetInterface< float> &FEATURE) const
    {
      // predict on an unnormalized feature vector
      const FeatureDataSet< float> result( PredictWithoutRescaling( m_RescaleInput->operator()( FEATURE)));
      return m_RescaleOutput->DeScale( result);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read SupportVectorMachineMultiOutput from std::istream
    //! @param ISTREAM input stream containing object
    //! @return modified input stream
    std::istream &SupportVectorMachineMultiOutput::Read( std::istream &ISTREAM)
    {
      // read member
      io::Serialize::Read( m_Beta, ISTREAM);
      io::Serialize::Read( m_Bias, ISTREAM);
      io::Serialize::Read( m_Kernel, ISTREAM);
      io::Serialize::Read( m_NumberSupportVectors, ISTREAM);
      io::Serialize::Read( m_RescaleInput, ISTREAM);
      io::Serialize::Read( m_RescaleOutput, ISTREAM);
      io::Serialize::Read( m_SupportVectors, ISTREAM);
      // end
      return ISTREAM;
    }

    //! @brief write SupportVectorMachineMultiOutput into std::ostream
    //! @param OSTREAM original output stream
    //! @param INDENT number of spaces for indent
    //! @return stream with appended data
    std::ostream &SupportVectorMachineMultiOutput::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member
      io::Serialize::Write( m_Beta, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Bias, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Kernel, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_NumberSupportVectors, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_RescaleInput, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_RescaleOutput, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_SupportVectors, OSTREAM, INDENT);
      // end
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer SupportVectorMachineMultiOutput::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
        (
          "see http://en.wikipedia.org/wiki/Support_vector_machine"
        );

      parameters.AddInitializer
      (
        "kernel",
        "kernel used to map pairs of features onto a hyperspace",
        io::Serialization::GetAgent( &m_Kernel),
        "RBF(gamma=0.5)"
      );

      return parameters;
    }

  } // namespace model
} // namespace bcl
