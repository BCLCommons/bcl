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
#include "model/bcl_model_support_vector_kernel_rbf.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"
#include "model/bcl_model_feature_reference.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {
    // add the interface to the set of known implementations
    const util::SiPtr< const util::ObjectInterface> SupportVectorKernelRBF::s_Instance
    (
      util::Enumerated< SupportVectorKernelBase>::AddInstance( new SupportVectorKernelRBF())
    );

    //! @brief constructor
    SupportVectorKernelRBF::SupportVectorKernelRBF() :
      m_Gamma( 1.0)
    {
    }

    //! @brief constructor
    SupportVectorKernelRBF::SupportVectorKernelRBF( const float &GAMMA) :
      m_Gamma( GAMMA)
    {
    }

    //! @brief virtual copy constructor
    SupportVectorKernelRBF *SupportVectorKernelRBF::Clone() const
    {
      return new SupportVectorKernelRBF( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &SupportVectorKernelRBF::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return a string to call the proper device kernel function for opencl
    //! @return a string that launches the opencl device function
    std::string SupportVectorKernelRBF::GetDeviceKernelFunctionCallString() const
    {
      return "RBFKernelDevice( feature, support_vector, length,   " + util::Format()( m_Gamma) + ")";
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &SupportVectorKernelRBF::GetAlias() const
    {
      static const std::string s_Name( "RBF");
      return s_Name;
    }

    //! @brief gets the scheme for this function class
    //! @return the scheme for this function class
    const std::string &SupportVectorKernelRBF::GetScheme() const
    {
      static const std::string s_scheme( "k( xi, xj)=e^( -gamma * Variance(xi, xj))");
      return s_scheme;
    }

    //! @brief calculate a kernel value for two vectors A and B
    //! @return a kernel value for two vectors A and B
    float SupportVectorKernelRBF::operator()
    (
      const FeatureReference< float> &VECTOR_A,
      const FeatureReference< float> &VECTOR_B
    ) const
    {
      float square_norm( 0.0);
      for
      (
        const float *vec_a( VECTOR_A.Begin()), *vec_b( VECTOR_B.Begin()), *vec_a_end( VECTOR_A.End());
        vec_a != vec_a_end;
        ++vec_a, ++vec_b
      )
      {
        square_norm += math::Sqr( *vec_a - *vec_b);
      }

      return std::exp( -m_Gamma * square_norm);
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer SupportVectorKernelRBF::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription( GetScheme());

      parameters.AddInitializer
      (
        "gamma",
        "coefficient in exponential of kernel, larger values attenuate the returned result",
        io::Serialization::GetAgent( &m_Gamma),
        "1.0"
      );

      return parameters;
    }

  } // namespace model
} // namespace bcl
