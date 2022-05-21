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

#ifndef BCL_MODEL_SUPPORT_VECTOR_MACHINE_H_
#define BCL_MODEL_SUPPORT_VECTOR_MACHINE_H_

// include the namespace header
#include "bcl_model.h"

// include other forward headers - sorted alphabetically
#include "math/bcl_math.fwd.hh"
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_model_feature_data_set.h"
#include "bcl_model_interface.h"
#include "bcl_model_support_vector_kernel_base.h"
#include "util/bcl_util_implementation.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SupportVectorMachine
    //! @brief is an support vector and is FunctionInterface derived
    //! @details It calculates for an argument vector (feature) of floats a prediction vector (result) with floats.
    //!
    //! @see @link example_model_support_vector_machine.cpp @endlink
    //! @author butkiem1
    //! @date 10/09/2009
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SupportVectorMachine :
      public Interface
    {
    public:

    //////////
    // data //
    //////////

      //! the default input range for SVM normalization range
      static const math::Range< float> s_DefaultInputRange;

    private:

    //////////
    // data //
    //////////

      //! @brief LaGrange Multipliers alpha
      storage::Vector< float> m_Alpha;

      //! @brief threshold of hyperplane
      float m_Bias;

      //! @brief Kernel that is used
      util::Implementation< SupportVectorKernelBase> m_Kernel;

      //! @brief number of BSV's
      size_t m_NumberBoundSupportVectors;

      //! @brief number of SV's
      size_t m_NumberSupportVectors;

      //! rescale input to normalization function input range
      util::ShPtr< RescaleFeatureDataSet> m_RescaleInput;

      //! rescale output from normalization function output range
      util::ShPtr< RescaleFeatureDataSet> m_RescaleOutput;

      //! @brief matrix where finally all support vectors are stored
      FeatureDataSet< float> m_SupportVectors;

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief standard constructor
      SupportVectorMachine();

      //! @brief standard constructor
      SupportVectorMachine
      (
        const util::Implementation< SupportVectorKernelBase> &KERNEL,
        const util::ShPtr< RescaleFeatureDataSet> &RESCALE_INPUT,
        const util::ShPtr< RescaleFeatureDataSet> &RESCALE_OUTPUT
      );

      //! @brief standard constructor
      SupportVectorMachine( const util::Implementation< SupportVectorKernelBase> &KERNEL);

      //! @brief Clone function
      //! @return pointer to current support vector machine
      SupportVectorMachine *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief get a short name for this class
      //! @return a short name for this class
      const std::string &GetAlias() const;

      //! @brief get the output feature size for this model
      //! @return the output feature size for this model
      size_t GetNumberOutputs() const;

    /////////////////
    // get and set //
    /////////////////

      //! @brief method for getting the support vectors
      //! @return util::ShPtr to DataSetInterface that contains the support vectors
      const FeatureDataSet< float> &GetSupportVectors() const
      {
        return m_SupportVectors;
      }

      //! @brief method for setting the support vectors
      void SetSupportVectors( const FeatureDataSet< float> &SUPPORTVECTORS)
      {
        m_SupportVectors = SUPPORTVECTORS;
      }

      //! @brief method for getting m_Kernel
      //! @return the used SupportVectorKernelBase
      const util::Implementation< SupportVectorKernelBase> &GetKernel() const
      {
        return m_Kernel;
      }

      //! @brief method for setting the m_Kernel
      //! @param KERNEL a specific SupportVectorKernelBase
      void SetKernel( const util::Implementation< SupportVectorKernelBase> &KERNEL)
      {
        m_Kernel = KERNEL;
      }

      //! @brief method for getting the m_Alpha
      storage::Vector< float> const &GetAlpha() const
      {
        return m_Alpha;
      }

      //! @param ALPHA a vector of LaGrange Multiplier
      void SetAlpha( const storage::Vector< float> &ALPHA)
      {
        m_Alpha = ALPHA;
      }

      //! @brief method for getting the m_Bias
      //! @return value of m_Bias
      const float GetBias() const
      {
        return m_Bias;
      }

      //! @brief method for setting the m_Bias
      //! @param BIAS
      void SetBias( const float BIAS)
      {
        m_Bias = BIAS;
      }

      //! @brief method for getting m_NumberBoundSupportVectors
      //! @return value of m_NumberBoundSupportVectors
      const size_t &GetNumberBoundSupportVectors() const
      {
        return m_NumberBoundSupportVectors;
      }

      //! @brief method for setting m_NumberSupportVectors
      //! @param NUMBER_BOUND_SUPPORT_VECTORS
      void SetNumberBoundSupportVectors( const size_t NUMBER_BOUND_SUPPORT_VECTORS)
      {
        m_NumberBoundSupportVectors = NUMBER_BOUND_SUPPORT_VECTORS;
      }

      //! @brief method for getting m_NumberBoundSupportVectors
      //! @return value of m_NumberSupportVectors
      const size_t &GetNumberSupportVectors() const
      {
        return m_NumberSupportVectors;
      }

      //! @brief method for setting m_NumberSupportVectors
      //! @param NUMBER_SUPPORT_VECTORS
      void SetNumberSupportVectors( const size_t NUMBER_SUPPORT_VECTORS)
      {
        m_NumberSupportVectors = NUMBER_SUPPORT_VECTORS;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief Set the scaling of a feature set according to the model
      //! @param FEATURES feature set of interest
      //! @note this allows for external classes that own a dataset to ensure that a new dataset is never created
      //!       when operator() is called
      void Rescale( FeatureDataSet< float> &FEATURE) const;

      //! @brief predict result with model using a NOT rescaled feature vector
      //! @param FEATURE not rescaled feature vector
      //! @return predcited result vector using a model
      FeatureDataSet< float> PredictWithoutRescaling( const FeatureDataSetInterface< float> &FEATURE) const;

      //! @brief predict result with model using a rescaled feature vector
      //! @param FEATURE normalized or rescaled feature vector
      //! @return predcited result vector using a model
      FeatureDataSet< float> operator()( const FeatureDataSetInterface< float> &FEATURE) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

      //! read NeuralNetwork from std::istream
      std::istream &Read( std::istream &ISTREAM);

      //! write NeuralNetwork into std::ostream
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class NeuralNetwork

  } // namespace model
} // namespace bcl

#endif // BCL_MODEL_SUPPORT_VECTOR_MACHINE_H_
