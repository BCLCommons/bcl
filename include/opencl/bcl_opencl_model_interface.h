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

#ifndef BCL_OPENCL_MODEL_INTERFACE_H_
#define BCL_OPENCL_MODEL_INTERFACE_H_

// include the namespace header
#include "bcl_opencl.h"

// include other forward headers - sorted alphabetically
#include "model/bcl_model.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_opencl_matrix.h"
#include "math/bcl_math_function_interface_serializable.h"
#include "model/bcl_model_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opencl
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ModelInterface
    //! @brief to unify all trainable model in the approximator framework that use gpus. This interface ensures the
    //!        possibility to predict a result for a normalized / rescaled feature vector as well as for a unnormalized
    //!        feature vector.
    //!
    //! @see @link example_opencl_model_interface.cpp @endlink
    //! @author loweew
    //! @date Apr 08, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ModelInterface :
      public model::Interface
    {

    public:

    ////////////////
    // operations //
    ////////////////

    ///////////////
    // operators //
    ///////////////

      //! @brief predict result with model without rescaling and using data already on the gpu
      //! @param FEATURE normalized or rescaled feature vector
      //! @return predicted results on the gpu
      virtual Matrix< float> operator()( const Matrix< float> &FEATURE) const = 0;

      //! @brief predict result with model without rescaling and using data already on the gpu
      //! @param FEATURE normalized or rescaled feature vector
      //! @return predicted results on the gpu
      virtual model::FeatureDataSet< float> operator()( const model::FeatureDataSetInterface< float> &FEATURE) const = 0;

    }; // class Interface

  } // namespace opencl
} // namespace bcl

#endif // BCL_OPENCL_MODEL_INTERFACE_H_
