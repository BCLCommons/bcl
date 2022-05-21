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

#ifndef BCL_MODEL_OBJECTIVE_FUNCTION_WRAPPER_H_
#define BCL_MODEL_OBJECTIVE_FUNCTION_WRAPPER_H_

// include the namespace header
#include "bcl_model.h"

// include other forward headers - sorted alphabetically
#include "linal/bcl_linal.fwd.hh"
#include "storage/bcl_storage.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_model_feature_data_set.h"
#include "bcl_model_objective_function_interface.h"
#include "descriptor/bcl_descriptor_dataset.h"
#include "math/bcl_math_binary_function_interface_serializable.h"
#include "math/bcl_math_function_interface_serializable.h"
#include "opti/bcl_opti_improvement_type.h"
#include "storage/bcl_storage_pair.h"
#include "storage/bcl_storage_vector.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ObjectiveFunctionWrapper
    //! @brief to unify all trainable objective functions for model::Interfaces in
    //!        the approximator framework. It enforces that a monitoring dataset can be set and retrieved used in
    //!        the derived objective function.
    //!
    //! @see @link example_model_objective_function_wrapper.cpp @endlink
    //! @author mendenjl
    //! @date May 03, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ObjectiveFunctionWrapper :
      public math::FunctionInterfaceSerializable< util::PtrInterface< Interface>, float>
    {

    private:

    //////////
    // data //
    //////////

      //! @brief data set that is used for testing the progress of the model::Interface
      mutable util::ShPtr< descriptor::Dataset> m_Data;

      //! @brief specify objective function implementation
      util::Implementation< ObjectiveFunctionInterface> m_Specialization;

    public:

      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      ObjectiveFunctionWrapper();

      //! @brief constructor with parameters
      //! @param SPECIALIZATION specify objective function implementation
      explicit ObjectiveFunctionWrapper
      (
        const util::Implementation< ObjectiveFunctionInterface> &SPECIALIZATION
      );

      //! @brief constructor with parameters
      //! @param DATA data set that is used for testing the progress of the model::Interface
      //! @param SPECIALIZATION specify objective function implementation
      ObjectiveFunctionWrapper
      (
        util::ShPtr< descriptor::Dataset> &DATA,
        const util::Implementation< ObjectiveFunctionInterface> &SPECIALIZATION
      );

      //! @brief Clone function
      //! @return pointer to new Interface
      ObjectiveFunctionWrapper *Clone() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief set data and rescaling for support vector model
      //! @param DATA monitoring dataset
      //! @param RESCALE_INPUT the rescale function
      void SetData
      (
        util::ShPtr< descriptor::Dataset> &DATA,
        const util::ShPtr< RescaleFeatureDataSet> &RESCALE_INPUT
      );

      //! @brief get training data for support vector model
      const util::ShPtr< descriptor::Dataset> &GetData() const;

      //! @brief determine what sign of the derivative of this objective function indicates improvement
      //! @return the sign of the derivative of this objective function indicates improvement
      opti::ImprovementType GetImprovementType() const;

      //! @brief get the overall goal of the objective function
      //! @return the goal of the objective function
      ObjectiveFunctionInterface::Goal GetGoalType() const;

      //! @brief get implementation used by the objective function
      const util::Implementation< ObjectiveFunctionInterface> &GetImplementation() const;

      //! @brief get the threshold, for classification type objectives
      //! @return the threshold, for classification type objectives
      float GetThreshold() const;

      //! @brief get the parity, for rank classification type objectives
      //! @return the parity, for rank classification type objectives
      //! @note this has meaning only for rank-classification objectives; true means the objective is most interested
      //!       in prediction of values higher than the threshold, false means below threshold
      bool GetRankingParity() const;

      //! @brief get the desired hit rate (e.g. fraction of positive predictions)
      //! @return the desired hit rate, for rank classification type objectives
      //! @note this has meaning only for rank-classification objectives; 0.01 means, for example, that only the top 1%
      //!       of values will be considered
      float GetDesiredHitRate() const;

      //! @brief test whether one objective function value is better than another
      //! @param NEW_RESULT, OLD_RESULT two objective function results
      //! @return true if NEW_RESULT is better than OLD_RESULT
      bool TestWhetherResultsImproved( const float &NEW_RESULT, const float &OLD_RESULT);

      //! @brief alter an output descaling function based on the objective function
      //! @param RESCALING current rescaling function
      //! @param PREDICTIONS actual model predictions
      //! @return the optimized rescaling function
      util::ShPtr< RescaleFeatureDataSet> OptimizeRescalingFunction
      (
        const util::ShPtr< RescaleFeatureDataSet> &RESCALE,
        const FeatureDataSet< float> &PREDICTIONS
      ) const;

    ///////////////
    // operators //
    ///////////////

      //! @brief compute error between predicted values from FEATURE and RESULT
      //! @param MODEL shptr on model::Interface that will be monitored
      //! @return error between predicted values from FEATURE and RESULT
      float operator()( const util::PtrInterface< Interface> &MODEL) const;

    //////////
    // read //
    //////////

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief compute a list of pairs with experimental and predicted values
      //! @param MODEL shptr on model::Interface that will be monitored
      //! @return feature data set of predicted results
      FeatureDataSet< float> Predict( const util::PtrInterface< Interface> &MODEL) const;

      //! @brief compute feature dataset with predicted values
      //! @param PREDICTIONS predictions for which to evaluate this measure
      //!        PREDICTIONS may be rescaled, as necessary, for this operation, but it will be returned
      //!        with the same scaling it was given during input
      //! @return predicted values from FEATURE and RESULT based on given model
      float Evaluate( FeatureDataSet< float> &PREDICTIONS) const;

    }; // class ObjectiveFunctionWrapper
  } // namespace model
} // namespace bcl

#endif // BCL_MODEL_OBJECTIVE_FUNCTION_WRAPPER_H_
