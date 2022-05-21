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

#ifndef BCL_MODEL_INTERFACE_H_
#define BCL_MODEL_INTERFACE_H_

// include the namespace header
#include "bcl_model.h"

// include other forward headers - sorted alphabetically
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_function_interface_serializable.h"
#include "util/bcl_util_enumerated.h"
#include "util/bcl_util_serializable_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Interface
    //! @brief to unify all trainable model in the approximator framework. This interface ensures the
    //!        possibility to predict a result for a normalized / rescaled feature vector as well as for a unnormalized
    //!        feature vector.
    //!
    //! @see @link example_model_interface.cpp @endlink
    //! @author butkiem1, loweew
    //! @date Feb 24, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Interface :
      public math::FunctionInterfaceSerializable< FeatureDataSetInterface< float>, FeatureDataSet< float> >
    {

    public:

    /////////////////
    // data access //
    /////////////////

      //! @brief virtual copy constructor
      virtual Interface *Clone() const = 0;

      //! @brief get the output feature size for this model
      //! @return the output feature size for this model
      virtual size_t GetNumberOutputs() const = 0;

    ////////////////
    // operations //
    ////////////////

      //! @brief predict result with model using a NOT rescaled feature vector
      //! @param FEATURE not rescaled feature vector
      //! @return predcited result vector using a model
      virtual FeatureDataSet< float> PredictWithoutRescaling( const FeatureDataSetInterface< float> &FEATURE) const = 0;

      //! @brief Set the scaling of a feature set according to the model
      //! @param FEATURES feature set of interest
      //! @note this allows for external classes that own a dataset to ensure that a new dataset is never created
      //!       when operator() is called
      virtual void Rescale( FeatureDataSet< float> &FEATURE) const = 0;

    ///////////////
    // operators //
    ///////////////

      //! @brief predict result with model using a rescaled feature vector
      //! @param FEATURE normalized or rescaled feature vector
      //! @return predcited result vector using a model
      virtual FeatureDataSet< float> operator()( const FeatureDataSetInterface< float> &FEATURE) const = 0;

    }; // class Interface

  } // namespace model
} // namespace bcl

#endif // BCL_MODEL_INTERFACE_H_
