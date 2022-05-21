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

#ifndef BCL_MODEL_SUPPORT_VECTOR_KERNEL_BASE_H_
#define BCL_MODEL_SUPPORT_VECTOR_KERNEL_BASE_H_

// include the namespace header
#include "bcl_model.h"

// include other forward headers - sorted alphabetically
#include "linal/bcl_linal.fwd.hh"
#include "storage/bcl_storage.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_enumerated.h"
#include "util/bcl_util_serializable_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SupportVectorKernelBase
    //! @brief Base Class for Kernel classes to provide a common structure among all kernels. This interface defines
    //!        which methods are required in a kernel class in order to fit into the kernel scheme.
    //!
    //! @remarks example unnecessary
    //! @author butkiem1
    //! @date 04.08.2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SupportVectorKernelBase :
      public util::SerializableInterface
    {

    public:

    /////////////////
    // data access //
    /////////////////

      //! @return new Pointer to a copy of the actual object behind the pointer
      virtual SupportVectorKernelBase *Clone() const = 0;

    ////////////////
    // operations //
    ////////////////

      //! @brief return a storage::Vector of calculated values for a given SmallMolecule object
      virtual float operator()
      (
        const FeatureReference< float> &VECTOR_A,
        const FeatureReference< float> &VECTOR_B
      ) const = 0;

      //! @brief compute kernel values for every pair input_vector_i with each other vector in data set
      //! @param TRAINING_DATA  feature vector data set without labels
      //! @param INPUT_VECTOR_I feature vector
      //! @return vector with all combinations towards INPUT_VECTOR_I
      linal::Vector< float> GetInputVectorIKernelMatrixMultiOutput
      (
        const FeatureDataSetInterface< float> &TRAINING_DATA,
        const size_t &INPUT_VECTOR_I
      ) const;

      //! @brief compute kernel values for every pair input_vector_i with each other vector in data set
      //! @param TRAINING_DATA  feature vector data set without labels
      //! @param INPUT_VECTOR_I feature vector
      //! @param SIGNS
      //! @param PROBLENGTH
      //! @return vector with all combinations towards INPUT_VECTOR_I
      linal::Vector< float> GetInputVectorIKernelMatrix
      (
        const FeatureDataSetInterface< float> &TRAINING_DATA,
        const size_t &INPUT_VECTOR_I,
        const storage::Vector< int> &SIGNS,
        const size_t &PROBLENGTH
      ) const;

      //! @brief return a string to call the proper device kernel function for opencl
      //! @return a string of the form "RBFKernelDevice( feature, support_vector, length, 2.2)"
      virtual std::string GetDeviceKernelFunctionCallString() const = 0;

    private:

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief compute kernel matrix for one particular feature with all feature values in training data
      //! @param TRAIN_INDEX_RESULT data structure which contains pointer to training data, current vector,
      //!        and result vector
      void ComputePartialInputVectorIKernelMatrix
      (
        storage::Triplet
        <
          util::SiPtr< const FeatureDataSetInterface< float> >,
          const FeatureReference< float>,
          linal::Vector< float>
        > &TRAIN_INDEX_RESULT,
        const size_t &START_INDEX,
        const size_t &END_INDEX
      ) const;

    }; // class SupportVectorKernelBase

  } // namespace model
} // namespace bcl

#endif // BCL_MODEL_SUPPORT_VECTOR_KERNEL_BASE_H_
