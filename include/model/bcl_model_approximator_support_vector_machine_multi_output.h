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

#ifndef BCL_MODEL_APPROXIMATOR_SUPPORT_VECTOR_MACHINE_MULTI_OUTPUT_H_
#define BCL_MODEL_APPROXIMATOR_SUPPORT_VECTOR_MACHINE_MULTI_OUTPUT_H_

// include the namespace header
#include "bcl_model.h"

// include other forward headers - sorted alphabetically
#include "math/bcl_math.fwd.hh"
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_model_approximator_base.h"
#include "bcl_model_feature_data_set.h"
#include "bcl_model_objective_function_wrapper.h"
#include "descriptor/bcl_descriptor_dataset.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ApproximatorSupportVectorMachineMultiOutput
    //! @brief provides the iteration function for training a
    //!        support vector machine with the sequential minimal optimization algorithm with multiple outputs
    //!
    //! @see @link example_model_approximator_support_vector_machine_multi_output.cpp @endlink
    //! @author krauseg, butkiem1
    //! @date Apr 16, 2014
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ApproximatorSupportVectorMachineMultiOutput :
      public ApproximatorBase
    {
    private:

    //////////
    // data //
    //////////

      //! @brief indicates upper boundary
      static const size_t s_Upper_Bound;

      //! @brief indicates lower boundary
      static const size_t s_Lower_Bound;

      //! @brief indicates no boundary
      static const size_t s_Free;

      //! @brief epsilon as constant Value
      static const float m_EPS_A; //TODO rename to s_EPS_A

      //! @brief penalty parameter C
      float m_CostParameterC;

      //! @brief Status of feature vectors
      storage::Vector< storage::Vector< size_t> > m_Status;

      //! @brief LaGrange multipliers alpha
      storage::Vector< linal::Vector< float> > m_Alpha;

      //! @brief Gradients
      storage::Vector< linal::Vector< float> > m_Gradient;

      //! @ brief Bias/Offset
      linal::Vector< float> m_Bias;

      //! rescale input to transfer function input range
      util::ShPtr< RescaleFeatureDataSet> m_RescaleInput;

      //! rescale output from transfer function output range
      util::ShPtr< RescaleFeatureDataSet> m_RescaleOutput;

      //! @brief dataset with training data
      util::ShPtr< SupportVectorMachineMultiOutput> m_Model;

      //! @brief pointer to tracker that keeps track of results
      util::ShPtr< ObjectiveFunctionWrapper> m_ObjectiveFunction;

      //! @brief < number of iterations until objective function evaluation
      size_t m_NumberIterations;

      //! @brief < number of current support vectors
      size_t m_NumberCurrentSupportVectors;

      //! @brief < optimization gap threshold to stop training process
      float m_OptimizationGapThreshold;

      //! @brief < current optimization gap
      float m_OptimizationGap;

      //! @brief last objective function value
      float m_LastObjectiveFunctionValue;

      //! @brief single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

       //! @brief number of nearest neighbors kappa of a result; used for building LocallyLinearTransformation coordinate system
      size_t m_Kappa;

      //! @brief LLT coordinate systems for all result vectors
      storage::Vector< linal::Matrix< float> > m_CoordinateSystems;

      //! @brief Transposed LLT coordinate systems for all result vectors
      storage::Vector< linal::Matrix< float> > m_CoordinateSystemsTransposed;

      //! @brief if alpha is smaller than this threshold it is set to zero for the final model
      float m_Epsilon;

      //! @brief kernel value of the corresponding feature vector with itself
      linal::Vector< float> m_SelfKernelValues;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      ApproximatorSupportVectorMachineMultiOutput();

      //! @brief Iterate for Sequential Minimal Optimization Learning Algorithm
      //! @param COST_PARAMETER_C penalty parameter c for svm regression training
      //! @param MODEL initial support vector model
      //! @param TRAINING_DATA training data set of choice
      //! @param NUMBER_ITERATIONS
      ApproximatorSupportVectorMachineMultiOutput
      (
        const float COST_PARAMETER_C,
        const util::ShPtr< SupportVectorMachineMultiOutput> &MODEL,
        util::ShPtr< descriptor::Dataset> &TRAINING_DATA,
        const size_t NUMBER_ITERATIONS
      );

      //! @brief copy constructor
      //! @return a new ApproximatorSupportVectorMachineMultiOutput copied from this instance
      ApproximatorSupportVectorMachineMultiOutput *Clone() const;

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
      void SetTrainingData( util::ShPtr< descriptor::Dataset> &DATA);

    ////////////////
    // operations //
    ////////////////

      //! @brief construct a model from the current iterate
      //! @return shptr to the new model interface
      util::ShPtr< Interface> GetCurrentModel() const;

      //! @brief method for getting the m_CostParameterC
      //! @return penalty or cost parameter C
      float GetCostParameterC() const
      {
        return m_CostParameterC;
      }

      //! @brief method for setting the m_CostParameterC
      //! @param COST_PARAMETER_C penalty or cost parameter C
      void SetCostParameterC( const float COST_PARAMETER_C)
      {
        m_CostParameterC = COST_PARAMETER_C;
      }

      //! @brief returns the current approximation
      //! @return current argument result pair
      const util::ShPtr< storage::Pair< util::ShPtr< Interface>, float> > GetCurrentApproximation() const;

      //! @brief conducts the next approximation step and stores the approximation
      void Next();

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read ApproximatorSupportVectorMachineMultiOutput from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write ApproximatorSupportVectorMachineMultiOutput into std::ostream
      //! @param OSTREAM output stream
      //! @param INDENT indentation
      //! @return ostream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    private:

      //! @brief initialize a default objective function constant as default objective function
      //! @return Wrapper with the objective Function
      static util::ShPtr< ObjectiveFunctionWrapper> GetDefaultObjectiveFunction();

      //! @brief evaluates whether the approximation can continue
      //! @return true, if the approximation can continue - otherwise false
      bool CanContinue() const;

      //! @brief finalize support vector model and determine final SVs, Betas and Bias
      //! @return optimization gap
      float FinalizeSupportVectorModel();

      //! @brief iterates one cycle and returns ShPtr to pair of resultant argument and corresponding score
      float IterationStep();

      //! @brief update status of data vectors according to alpha value of specific vector
      //! @param FEATURE_VECTOR_I
      void UpdateAlphaStatus( const size_t &FEATURE_VECTOR_I);

      //! @brief method checks whether a feature_vector i reached the upper boundary
      //! @param FEATURE_VECTOR_I a position of a feature vector in m_Status
      //! @param ELEMENT position of the element - corresponds to the output column
      //! @return const bool that indicates whether FEATURE_VECTOR_I reached upper boundary
      bool IsUpperBound( const size_t &FEATURE_VECTOR_I, const size_t &ELEMENT) const;

      //! @brief checks whether a feature_vector i reached the lower boundary
      //! @param FEATURE_VECTOR_I a position of a feature vector in m_Status
      //! @param ELEMENT position of the element - corresponds to the output column
      //! @return const bool that indicates whether FEATURE_VECTOR_I reached lower boundary
      bool IsLowerBound( const size_t &FEATURE_VECTOR_I, const size_t &ELEMENT) const;

      //! @brief checks whether a feature_vector i is not bound and between the boundaries
      //! @param FEATURE_VECTOR_I a position of a feature vector in m_Status
      //! @param ELEMENT position of the element - corresponds to the output column
      //! @return const bool that indicates whether FEATURE_VECTOR_I reached no boundary
      bool IsFree( const size_t &FEATURE_VECTOR_I, const size_t &ELEMENT) const;

      //! @brief Initialize member vectors and variables to prepare SVR training
      void InitializeMemberVectorsForTraining();

      //! @brief determine the bias vector for Support Vector Regression
      //! @return bias vector
      linal::Vector< float> CalculateBias() const;

      //! @brief heuristic approach to find a feature_vector i and j
      //! @brief depending on m_Gradients for examination of SMO Classification
      //! @param FIRST_VECTOR_INDEX index of first picked vector i
      //! @param SECOND_VECTOR_INDEX index of second picked vector j
      void DetermineFeatureVectorCombination
      (
        size_t &FIRST_VECTOR_INDEX, //TODO i, j
        size_t &SECOND_VECTOR_INDEX
      );

      //! @brief solve sub problem of quadratic problem according two feature vectors
      //! @param FEATURE_VECTOR_I index of feature vector in feature matrix
      //! @param FEATURE_VECTOR_J index of feature vector in feature matrix
      //! @param APPLY_SOLUTION indicates whether the solution is to be applied to m_Alpha, m_Status, m_Gradients
      //! @return absolute value of delta beta
      float SolveQuadraticProblemSubProblem
      (
        const size_t &FEATURE_VECTOR_I,
        const size_t &FEATURE_VECTOR_J, //TODO size_t
        const bool &APPLY_SOLUTION
      );

      //! @brief computes for each output vector a local coordinate system used in the multi output method
      void ComputeCoordinateSystems();

    }; // class ApproximatorSupportVectorMachineMultiOutput

  } // namespace model
} // namespace bcl

#endif // BCL_MODEL_APPROXIMATOR_SUPPORT_VECTOR_MACHINE_MULTI_OUTPUT_H_
