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

#ifndef BCL_MODEL_SUPPORT_VECTOR_KERNEL_POLYNOMIAL_H_
#define BCL_MODEL_SUPPORT_VECTOR_KERNEL_POLYNOMIAL_H_

// include the namespace header
#include "bcl_model.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_model_support_vector_kernel_base.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SupportVectorKernelPolynomial
    //! @brief class for the POLYNOMIAL support vector kernel. The implementation and helper functions of
    //!        that specific kernel are realized here.
    //!        k(xi,xj) = (xi * xj + h)^d // * is the dot product; h is the homogeneity; d is the exponent
    //! in order to get a linear kernel, use an exponent of 1, and homogeneity of 1
    //!
    //!
    //! @see @link example_model_support_vector_kernel_polynomial.cpp @endlink
    //! @author butkiem1
    //! @date Aug 1, 2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SupportVectorKernelPolynomial :
      public SupportVectorKernelBase
    {

    private:

      //! exponent of polynomial kernel
      size_t m_Exponent;

      //! homogeneity
      bool m_Homogeneous;

      //! whether this class has been registered with the known implementations
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      SupportVectorKernelPolynomial() :
        m_Exponent( size_t( 1)),
        m_Homogeneous( true)
      {
      }

      //! @brief virtual copy constructor
      SupportVectorKernelPolynomial *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief gets the scheme for this function class
      //! @return the scheme for this function class
      const std::string &GetScheme() const;

      //! @brief return a string to call the proper device kernel function for opencl
      //! @return a string that launches the opencl device function
      std::string GetDeviceKernelFunctionCallString() const;

      //! @brief get a short name for this class
      //! @return a short name for this class
      const std::string &GetAlias() const;

      //! @brief calculate a kernel value for two vectors A and B
      //! @return a kernel value for two vectors A and B
      float operator()
      (
        const FeatureReference< float> &VECTOR_A,
        const FeatureReference< float> &VECTOR_B
      ) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      virtual io::Serializer GetSerializer() const;

    }; // class SupportVectorKernelPolynomial

  } // namespace model
} // namespace bcl

#endif // BCL_MODEL_SUPPORT_VECTOR_KERNEL_POLYNOMIAL_H_
