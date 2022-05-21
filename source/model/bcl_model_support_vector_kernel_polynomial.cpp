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
#include "model/bcl_model_support_vector_kernel_polynomial.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"
#include "linal/bcl_linal_vector_operations.h"
#include "model/bcl_model_feature_reference.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {
    // add the interface to the set of known implementations
    const util::SiPtr< const util::ObjectInterface> SupportVectorKernelPolynomial::s_Instance
    (
      util::Enumerated< SupportVectorKernelBase>::AddInstance( new SupportVectorKernelPolynomial())
    );

    //! @brief clone
    //! @return pointer to cloned object
    SupportVectorKernelPolynomial *SupportVectorKernelPolynomial::Clone() const
    {
      return new SupportVectorKernelPolynomial( *this);
    }

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &SupportVectorKernelPolynomial::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return a string to call the proper device kernel function for opencl
    //! @return a string that launches the opencl device function
    std::string SupportVectorKernelPolynomial::GetDeviceKernelFunctionCallString() const
    {
      return "PolynomialKernelDevice( feature, support_vector, length,   " + util::Format()( m_Exponent) + ", " + util::Format()( m_Homogeneous) + ")";
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &SupportVectorKernelPolynomial::GetAlias() const
    {
      static const std::string s_name( "Polynomial");
      return s_name;
    }

    //! @brief gets the scheme for this function class
    //! @return the scheme for this function class
    const std::string &SupportVectorKernelPolynomial::GetScheme() const
    {
      static const std::string s_desciption( "k(xi, xj) = (xi * xj + homogeneity) ^ exponent");
      return s_desciption;
    }

    //! @brief calculate the kernel value for the linear kernel according to two input vectors
    float SupportVectorKernelPolynomial::operator()
    (
      const FeatureReference< float> &VECTOR_A,
      const FeatureReference< float> &VECTOR_B
    ) const
    {
      // scalar product  of two vectors
      return math::Pow< float>( linal::ScalarProduct( VECTOR_A, VECTOR_B) + float( m_Homogeneous), float( m_Exponent));
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer SupportVectorKernelPolynomial::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription( GetScheme());

      parameters.AddInitializer
      (
        "exponent",
        "exponent of the function; allows to model feature conjunctions",
        io::Serialization::GetAgentWithRange( &m_Exponent, 0, 10),
        "1"
      );
      parameters.AddInitializer
      (
        "homogeneity",
        "homogeneity of the polynomial - 0 for homogenous, 1 for inhomogenous",
        io::Serialization::GetAgent( &m_Homogeneous),
        "1"
      );

      return parameters;
    }

  } // namespace model
} // namespace bcl
