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

#ifndef BCL_MODEL_APPROXIMATOR_BASE_H_
#define BCL_MODEL_APPROXIMATOR_BASE_H_

// include the namespace header
#include "bcl_model.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_model_has_labels_base.h"
#include "bcl_model_interface.h"
#include "bcl_model_objective_function_wrapper.h"
#include "bcl_model_rescale_feature_data_set.h"
#include "descriptor/bcl_descriptor_dataset.h"
#include "opti/bcl_opti_approximator_modular_base.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ApproximatorBase
    //! @brief This class provides basic functionality for a approximation scheme in the model context that allows
    //! specialization for different problems.
    //!
    //! @remarks example unnecessary
    //! @author fischea
    //! @date Jan 5, 2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ApproximatorBase :
      public opti::ApproximatorModularBase< util::ShPtr< Interface>, float>,
      public util::SerializableInterface,
      virtual public HasLabelsBase
    {

    //////////
    // data //
    //////////

    protected:

      //! flag whether training was continued from a serialized file
      bool m_TrainingContinued;

      //! shared pointer to the training data
      util::ShPtr< descriptor::Dataset> m_TrainingData;

      //! shared pointer to the objective function used to evaluate a monitoring data set
      util::ShPtr< ObjectiveFunctionWrapper> m_ObjectiveFunction;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

    public:

      //! @brief default constructor
      //! @detail sets the improvement type of the tracker to SmallerIsBetter
      ApproximatorBase();

      //! @brief construct from objective function
      //! @detail sets the improvement type of the tracker to SmallerIsBetter
      //! @param OBJECTIVE_FUNCTION objective function used to evaluate a monitoring data set
      ApproximatorBase( const util::ShPtr< ObjectiveFunctionWrapper> &OBJECTIVE_FUNCTION);

      //! @brief Clone function
      //! @return pointer to a new ApproximatorBase
      virtual ApproximatorBase *Clone() const = 0;

    /////////////////
    // data access //
    /////////////////

      //! @brief construct a model from the current iterate
      //! @return shptr to the new model interface
      virtual util::ShPtr< Interface> GetCurrentModel() const = 0;

      //! @brief returns a shared pointer to the training data set
      //! @return shared pointer to the data set interface with training data
      const util::ShPtr< descriptor::Dataset> &GetTrainingData() const;

      //! @brief sets the training data set
      //! @param DATA training data set to be set
      virtual void SetTrainingData( util::ShPtr< descriptor::Dataset> &DATA) = 0;

      //! @brief get objective function to evaluate a monitoring data set
      //! @return data set interface with training data
      const util::ShPtr< ObjectiveFunctionWrapper> &GetObjectiveFunction() const;

      //! @brief sets the objective function used to evaluate a monitoring data set
      //! @param OBJECTIVE_FUNCTION objective function used to evaluate a monitoring data set
      virtual void SetObjectiveFunction( const util::ShPtr< ObjectiveFunctionWrapper> &OBJECTIVE_FUNCTION);

      //! @brief returns a shared pointer to the rescale function to transfer data sets into a scaled feature space
      //! @return shared pointer to rescale function to rescale into normalized feature space
      util::ShPtr< RescaleFeatureDataSet> GetRescaleFeatureDataSet() const;

      //! @brief returns a shared pointer to the rescale function to transfer data sets into a scaled result space
      //! @return shared pointer to the rescale function to rescale into normalized result space
      util::ShPtr< RescaleFeatureDataSet> GetRescaleResultDataSet() const;

      //! @brief returns the flag for continuation of training
      //! @return flag for continuation of training
      bool IsTrainingContinued() const;

      //! @brief sets the flag for continuation of training
      //! @param flag for continuation of training
      void SetTrainingContinued( const bool &CONTINUED_TRAINING);

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief Set the members of this property from the given LABEL
      //! @param LABEL the label to parse
      //! @param ERR_STREAM stream to write errors out to
      bool ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERR_STREAM);

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      virtual std::istream &Read( std::istream &ISTREAM) = 0;

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      virtual std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const = 0;

    }; // class ApproximatorBase

  } // namespace model
} // namespace bcl

#endif // BCL_MODEL_APPROXIMATOR_BASE_H_
