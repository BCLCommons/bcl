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
#include "model/bcl_model_support_vector_machine.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {

  //////////
  // data //
  //////////

    //! the default input range for SVM normalization range
    const math::Range< float> SupportVectorMachine::s_DefaultInputRange( 0, 1);

    // add the interface to the set of known implementations
    const util::SiPtr< const util::ObjectInterface> SupportVectorMachine::s_Instance
    (
      GetObjectInstances().AddInstance( new SupportVectorMachine())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief standard constructor
    SupportVectorMachine::SupportVectorMachine() :
      m_Alpha(),
      m_Bias(),
      m_Kernel(),
      m_NumberBoundSupportVectors(),
      m_NumberSupportVectors(),
      m_RescaleInput(),
      m_RescaleOutput(),
      m_SupportVectors()
    {
    }

    //! @brief standard constructor
    //! @param KERNEL
    //! @param RESCALE_INPUT
    //! @param RESCALE_OUTPUT
    SupportVectorMachine::SupportVectorMachine
    (
      const util::Implementation< SupportVectorKernelBase> &KERNEL,
      const util::ShPtr< RescaleFeatureDataSet> &RESCALE_INPUT,
      const util::ShPtr< RescaleFeatureDataSet> &RESCALE_OUTPUT
    ) :
      m_Alpha( storage::Vector< float>()),
      m_Bias( float( 0)),
      m_Kernel( KERNEL),
      m_NumberBoundSupportVectors( size_t( 0)),
      m_NumberSupportVectors( size_t( 0)),
      m_RescaleInput( RESCALE_INPUT),
      m_RescaleOutput( RESCALE_OUTPUT),
      m_SupportVectors( FeatureDataSet< float>( linal::MatrixConstReference< float>()))
    {
    }

    //! @brief standard constructor
    //! @param KERNEL
    SupportVectorMachine::SupportVectorMachine
    (
      const util::Implementation< SupportVectorKernelBase> &KERNEL
    ) :
      m_Alpha( storage::Vector< float>()),
      m_Bias( float( 0)),
      m_Kernel( KERNEL),
      m_NumberBoundSupportVectors( size_t( 0)),
      m_NumberSupportVectors( size_t( 0)),
      m_RescaleInput( new RescaleFeatureDataSet()),
      m_RescaleOutput( new RescaleFeatureDataSet()),
      m_SupportVectors( FeatureDataSet< float>())
    {
    }

    //! @brief virtual copy constructor
    SupportVectorMachine *SupportVectorMachine::Clone() const
    {
      return new SupportVectorMachine( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &SupportVectorMachine::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &SupportVectorMachine::GetAlias() const
    {
      static const std::string s_Name( "SupportVectorMachine");
      return s_Name;
    }

    //! @brief get the output feature size for this model
    //! @return the output feature size for this model
    size_t SupportVectorMachine::GetNumberOutputs() const
    {
      return size_t( 1);
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
    void SupportVectorMachine::Rescale( FeatureDataSet< float> &FEATURE) const
    {
      if( !FEATURE.IsRescaled() || *FEATURE.GetScaling() != *m_RescaleInput)
      {
        FEATURE.DeScale();
        FEATURE.Rescale( *m_RescaleInput);
      }
    }

    //! @brief predict with normalized data
    FeatureDataSet< float> SupportVectorMachine::PredictWithoutRescaling
    (
      const FeatureDataSetInterface< float> &FEATURE
    ) const
    {
      const size_t num_points( FEATURE.GetNumberFeatures());
      storage::Vector< float> pred;
      pred.AllocateMemory( num_points);

      for( size_t i( 0); i < num_points; ++i)
      {
        // initialize classification result with bias
        linal::Vector< float> classification( 1, -m_Bias);

        // progress counter for accessing data set with support vectors
        size_t progress( 0);

        //iterator over all LaGrange multiplier alpha
        for
        (
          storage::Vector< float>::const_iterator iter_alpha_begin( GetAlpha().Begin()),
          iter_alpha_end( GetAlpha().End());
          iter_alpha_begin != iter_alpha_end;
          ++iter_alpha_begin, ++progress
        )
        {
          // calculate classification value of alpha and kernel value of unknown vector and a support vector
          classification( 0) +=
            *iter_alpha_begin * m_Kernel->operator()( FEATURE( i), m_SupportVectors( progress));
        }
        pred.PushBack( classification( 0));
      }

      // classification result
      return FeatureDataSet< float>
      (
        linal::Matrix< float>( num_points, 1, linal::Vector< float>( pred.Begin(), pred.End()).Begin()),
        *m_RescaleOutput
      );
    }

    //! predict with not normalized data vector
    FeatureDataSet< float> SupportVectorMachine::operator()( const FeatureDataSetInterface< float> &FEATURE) const
    {
      // handle the case where rescaling is necessary
      if( !FEATURE.IsRescaled() || *FEATURE.GetScaling() != *m_RescaleInput)
      {
        FeatureDataSet< float> feature( FEATURE);
        feature.Rescale( *m_RescaleInput);
        return PredictWithoutRescaling( feature).DeScale();
      }

      // data is already rescaled
      return PredictWithoutRescaling( FEATURE).DeScale();
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! read NeuralNetwork from std::istream
    std::istream &SupportVectorMachine::Read( std::istream &ISTREAM)
    {
      // read member
      io::Serialize::Read( m_Alpha, ISTREAM);
      io::Serialize::Read( m_Bias, ISTREAM);
      io::Serialize::Read( m_Kernel, ISTREAM);
      io::Serialize::Read( m_NumberBoundSupportVectors, ISTREAM);
      io::Serialize::Read( m_NumberSupportVectors, ISTREAM);
      io::Serialize::Read( m_RescaleInput, ISTREAM);
      io::Serialize::Read( m_RescaleOutput, ISTREAM);
      io::Serialize::Read( m_SupportVectors, ISTREAM);

      // end
      return ISTREAM;
    }

    //! write NeuralNetwork into std::ostream
    std::ostream &SupportVectorMachine::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member
      io::Serialize::Write( m_Alpha, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Bias, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Kernel, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_NumberBoundSupportVectors, OSTREAM, INDENT) << '\n';
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
    io::Serializer SupportVectorMachine::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
        (
          "see http://en.wikipedia.org/wiki/Support_vector_machine"
        );

      parameters.AddInitializer
      (
        "kernel",
        "kernel used to map pairs of features onto a hyperplane",
        io::Serialization::GetAgent( &m_Kernel),
        "RBF(gamma=0.5)"
      );

      return parameters;
    }

  } // namespace model
} // namespace bcl
