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

#ifndef BCL_MODEL_APPROXIMATOR_ALIGN_CUTOFF_H_
#define BCL_MODEL_APPROXIMATOR_ALIGN_CUTOFF_H_

// include the namespace header
#include "bcl_model.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_model_approximator_base.h"
#include "bcl_model_feature_data_set.h"
#include "bcl_model_neural_network_update_weights_interface.h"
#include "bcl_model_objective_function_wrapper.h"
#include "bcl_model_transfer_function_interface.h"
#include "linal/bcl_linal_vector.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ApproximatorAlignCutoff
    //! @brief Wrapper that maps values given by the internal iterate back onto the original threshold, so as to optimize precision
    //! Adjusts the values such that the model's output can be used with the input cutoff for rank-classification purposes
    //! Results of different machine learning methods or cross validation runs should not normally be averaged, for the
    //! purposes of making a jury prediction or presenting the results of a ranked-classification model, without using
    //! this wrapper.   This iterate applies a linear function to map the values given by the internal iterate.  The
    //! function is created so as to maximize precision of predictions, while maintaining a similar hit rate to the training data
    //!
    //! @see @link example_model_approximator_align_cutoff.cpp @endlink
    //! @author butkiem1, mendenjl
    //! @date Jul 07, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ApproximatorAlignCutoff :
      public ApproximatorBase
    {
    private:

    //////////
    // data //
    //////////

      //! internal iterate used for training
      util::Implementation< ApproximatorBase> m_InternalIterate;

      //! experimental cutoff
      float m_DescaledExperimentalCutoff;

      //! rescaled experimental cutoff
      float m_ExperimentalCutoff;

      //! cutoff specific to the determined classification model
      float m_ModelCutoff;

      //! average value of those above the cutoff
      float m_ExperimentalAverageAboveCutoff;

      //! average value of those below the cutoff
      float m_ExperimentalAverageBelowCutoff;

      //! average value of those above the cutoff
      float m_ModelAverageAboveCutoff;

      //! average value of those below the cutoff
      float m_ModelAverageBelowCutoff;

      //! parity = 0 if smaller values of the output value are more interesting, 1 if larger values are
      bool m_PositivesAboveThreshold;

      //! number of active compounds in training dataset
      size_t m_NumberActives;

      //! last internally generated model; used to prevent recalculating the threshold when the model didn't change
      util::ShPtr< Interface> m_LastInternallyGeneratedModel;
      util::ShPtr< Interface> m_LastModel;
      float m_LastObjectiveFunctionValue;

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

      //! @brief get the numerical tolerance used in this class
      //! @return the numerical tolerance used in this class
      static const float &GetNumericalTolerance();

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      ApproximatorAlignCutoff();

      //! @brief copy constructor
      //! @return a new ApproximatorAlignCutoff copied from this instance
      ApproximatorAlignCutoff *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief get a short name for this class
      //! @return a short name for this class
      const std::string &GetAlias() const;

      //! @brief set training data set for a specific iterate in approximater framework
      //! @param DATA training data set
      void SetTrainingData
      (
        util::ShPtr< descriptor::Dataset> &DATA
      );

      //! @brief set objective function to evaluate a monitoring dataset
      //! @param OBJ objective function of interest
      void SetObjectiveFunction( const util::ShPtr< ObjectiveFunctionWrapper> &OBJ);

    ////////////////
    // operations //
    ////////////////

      //! @brief construct a model from the current iterate
      //! @return shptr to the new model interface
      util::ShPtr< Interface> GetCurrentModel() const;

      //! @brief returns the current approximation
      //! @return current argument result pair
      const util::ShPtr< storage::Pair< util::ShPtr< Interface>, float> > GetCurrentApproximation() const;

      //! @brief conducts the next approximation step and stores the approximation
      void Next();

    protected:

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

      //! @brief Set the members of this property from the given LABEL
      //! @param LABEL the label to parse
      //! @param ERROR_STREAM the stream to write errors to
      bool ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM);

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    private:

      //! @brief evaluates whether the approximation can continue
      //! @return true, if the approximation can continue - otherwise false
      bool CanContinue() const;

      //! @brief determine the adjusted cutoff that maximizes the precision on the classification model
      void DetermineAdjustedCutoff();

    }; // class IterateResilientProagation

  } // namespace model
} // namespace bcl

#endif // BCL_MODEL_APPROXIMATOR_ALIGN_CUTOFF_H_
