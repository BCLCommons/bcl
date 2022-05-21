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

#ifndef BCL_OPENCL_APPROXIMATOR_SEQUENTIAL_MINIMIAL_OPTIMIZATION_H_
#define BCL_OPENCL_APPROXIMATOR_SEQUENTIAL_MINIMIAL_OPTIMIZATION_H_

// include the namespace header
#include "bcl_opencl.h"

// include other forward headers - sorted alphabetically
#include "math/bcl_math.fwd.hh"
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_opencl_command_queue.h"
#include "bcl_opencl_kernel_source_interface.h"
#include "bcl_opencl_matrix.h"
#include "bcl_opencl_support_vector_machine.h"
#include "bcl_opencl_vector.h"
#include "descriptor/bcl_descriptor_dataset.h"
#include "model/bcl_model_approximator_base.h"
#include "model/bcl_model_feature_data_set.h"
#include "model/bcl_model_objective_function_wrapper.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opencl
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ApproximatorSequentialMinimialOptimization
    //! @brief provides the iteration function for training a
    //!        support vector machine with the sequential minimal optimization algorithm in opencl - optimized for gpu
    //!
    //! @see @link example_opencl_approximator_sequential_minimial_optimization.cpp @endlink
    //! @author loweew
    //! @date Apr 14, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ApproximatorSequentialMinimialOptimization :
      public model::ApproximatorBase
    {
    private:

    //////////
    // data //
    //////////

      //! @brief indicates upper boundary
      static const cl_uint s_Upper_Bound;

      //! @brief indicates lower boundary
      static const cl_uint s_Lower_Bound;

      //! @brief indicates no boundary
      static const cl_uint s_Free;

      //! @brief epsilon as constant Value
      static const float m_EPS_A;

      //!< @brief parameter p (bias calculation)
      static const float m_P;

      //! @brief penalty parameter C
      float m_CostParameterC;

      //! @brief Status of feature vectors
      linal::Vector< cl_uint> m_Status;

      //! status on device
      Vector< cl_uint> m_StatusOnDevice;

      //! @brief LaGrange multipliers alpha
      linal::Vector< float> m_Alpha;

      //! alphas on device
      Vector< float> m_AlphaOnDevice;

      //! @brief Gradients
      linal::Vector< float> m_Gradient;

      //! gradients on device
      Vector< float> m_GradientOnDevice;

      //! @brief backup Gradients
      linal::Vector< float> m_GradientBar;

      //! backup of gradients?
      Vector< float> m_GradientBarOnDevice;

      //! @brief Bias / threshold
      linal::Vector< float> m_Bias;

      //! bias / threshold
      Vector< float> m_BiasOnDevice;

      //! kernel vectors
      Vector< float> m_Q_i;
      Vector< float> m_Q_j;
      Vector< float> m_Q_d;

      //! @brief vector of feature vector signs
      linal::Vector< cl_int> m_Signs;

      //! feature vector signs
      Vector< cl_int> m_SignsOnDevice;

      //! @brief vector of computed labels
      linal::Vector< float> m_Labels;

      //! vector of computed labels
      Vector< float> m_LabelsOnDevice;

      //! @brief determine working set of feature vectors
      size_t m_ActiveSize;

      //!< @brief length of the optimization problem
      size_t m_ProbLength;

      //! queue
      CommandQueue m_Queue;

      //! optimization gap threshold
      float m_OptimizationGapThreshold;

      //! optimization gap
      float m_OptimizationGap;

      //! @brief dataset with training data
      util::ShPtr< SupportVectorMachine> m_Model;

      //! training features on device
      Matrix< float> m_TrainingFeaturesOnDevice;

      //! training results on device
      Matrix< float> m_TrainingResultsOnDevice;

      //!< number of iterations until objective function evaluation
      size_t m_NumberIterations;

      //!< number of current support vectors
      size_t m_NumberCurrentSupportVectors;

      //! program
      cl::Program m_Program;

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

      //! block size
      static const cl_uint s_BlockSize = 16;

      float m_Gamma;

    public:

    //////////
    // data //
    //////////

      static const char *s_CLCompilerOptions;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      ApproximatorSequentialMinimialOptimization();

      //! @brief Iterate for Sequential Minimal Optimization Learning Algorithm
      //! @param COST_PARAMETER_C penalty parameter c for svm regression training
      //! @param OPENCL initial support vector model
      //! @param TRAINING_DATA training data set of choice
      //! @param NUMBER_ITERATIONS
      ApproximatorSequentialMinimialOptimization
      (
        const float COST_PARAMETER_C,
        const util::ShPtr< SupportVectorMachine> &MODEL,
        util::ShPtr< descriptor::Dataset> &TRAINING_DATA,
        const size_t NUMBER_ITERATIONS,
        const CommandQueue &QUEUE
      );

      //! @brief copy constructor
      //! @return a new ApproximatorSequentialMinimialOptimization copied from this instance
      ApproximatorSequentialMinimialOptimization *Clone() const;

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

      //! @brief set gamma
      //! @param GAMMA the gamma value
      void SetGamma( const float &GAMMA)
      {
        m_Gamma = GAMMA;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief finalize support vector model and determine final SVs, Alphas and Bias
      void FinalizeSupportVectorModel();

      //! @brief construct a model from the current iterate
      //! @return shptr to the new model interface
      util::ShPtr< model::Interface> GetCurrentModel() const;

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
      const util::ShPtr< storage::Pair< util::ShPtr< model::Interface>, float> > GetCurrentApproximation() const;

      //! @brief conducts the next approximation step and stores the approximation
      void Next();

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read ApproximatorSequentialMinimialOptimization from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write ApproximatorSequentialMinimialOptimization into std::ostream
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

      //! @brief evaluates whether the approximation can continue
      //! @return true, if the approximation can continue - otherwise false
      bool CanContinue() const;

      //! @brief iterates one cycle and returns ShPtr to pair of resultant argument and corresponding score
      float IterationStep();

      //! @brief method checks whether a LaGrange Multiplier reached a certain boundary or not
      //! @param ALPHA a LaGrange Multiplier
      //! @return const int that indicates whether ALPHA reached a certain boundary or not
      int AlphaToStatus( const float &ALPHA) const;

      //! @brief update status of data vectors according to alpha value of specific vector
      //! @param FEATURE_VECTOR_I
      void UpdateAlphaStatus( const int &FEATURE_VECTOR_I);

      //! @brief method checks whether a feature_vector i reached the upper boundary
      //! @param FEATURE_VECTOR_I a position of a feature vector in m_Status
      //! @return const bool that indicates whether FEATURE_VECTOR_I reached upper boundary
      bool IsUpperBound( const size_t &FEATURE_VECTOR_I) const;

      //! @brief checks whether a feature_vector i reached the lower boundary
      //! @param FEATURE_VECTOR_I a position of a feature vector in m_Status
      //! @return const bool that indicates whether FEATURE_VECTOR_I reached lower boundary
      bool IsLowerBound( const size_t &FEATURE_VECTOR_I) const;

      //! @brief checks whether a feature_vector i is not bound and between the boundaries
      //! @param FEATURE_VECTOR_I a position of a feature vector in m_Status
      //! @return const bool that indicates whether FEATURE_VECTOR_I reached no boundary
      bool IsFree( const size_t &FEATURE_VECTOR_I) const;

      //! @brief Initialize member vectors and variables to prepare SVR training
      //! @param TRAINING_DATA feature vector data set with labels
      //! @param SVR_OPENCL model to be trained
      void InitializeMemberVectorsForTraining();

      //! @brief determine the bias value for Support Vector Regression
      //! @return bias value
      float CalculateBias();

      //! @brief
      //! @param TRAINING_DATA
      //! @param SVR_OPENCL
      void ReconstructGradient();

      //! @brief heuristic approach to find a feature_vector i and j
      //! @brief depending on m_Gradients for examination of SMO Classification
      //! @param TRAINING_DATA data set of feature vectors inclusive labels
      //! @param SVR_OPENCL
      //! @param FIRST_VECTOR_INDEX
      //! @param SECOND_VECTOR_INDEX
      //! @return optimization gap
      float DetermineFeatureVectorCombination
      (
        size_t &FIRST_VECTOR_INDEX,
        size_t &SECOND_VECTOR_INDEX
      );

      //! @brief solve sub problem of quadratic problem according two feature vectors
      //! @param TRAINING_DATA training data set
      //! @param SVR_OPENCL support vector machine model
      //! @param FEATURE_VECTOR_I index of feature vector in feature matrix
      //! @param FEATURE_VECTOR_J index of feature vector in feature matrix
      void SolveQuadraticProblemSubProblem
      (
        const int &FEATURE_VECTOR_I,
        const int &FEATURE_VECTOR_J
      );

    //////////////////////
    // opencl functions //
    //////////////////////

      void oclFindMaxGradient
      (
        float &MAX_GRAD,
        cl_int &MAX_IND
      );

      void oclGetInputIKernelMatrixResultingVector
      (
        const cl_uint &VECTOR_ID,
        const float &GAMMA,
        Vector< float> &OUTPUT,
        Vector< float> &KERNEL_VECTOR
      );

      void oclGetGminGmax2
      (
        Vector< float> &Q_I,
        Vector< float> &Q_D,
        const cl_int &INDEX_I,
        float &MAX_GRADIENT,
        Vector< float> &TMP_GMAX2_OUTPUT,
        Vector< float> &TMP_OBJ_DIFF_MIN_OUTPUT,
        Vector< cl_uint> &TMP_GMIN_INDEX_OUTPUT,
        float &FINAL_GMAX2,
        cl_int       &FINAL_GMIN_IND
      );

      void oclCalculateBias
      (
        Vector< float> &UPPER_OUTPUT,
        Vector< float> &LOWER_OUTPUT,
        Vector< cl_uint>      &NR_FREE_OUTPUT,
        Vector< float> &SUM_FREE_OUTPUT,
        float &FINAL_UPPER_OUTPUT,
        float &FINAL_LOWER_OUTPUT,
        cl_uint      &FINAL_NR_FREE_OUTPUT,
        float &FINAL_SUM_FREE_OUTPUT
      );

      void oclUpdateAllAlphaStatus
      (
      );

      void oclUpdateGradient
      (
        const float &DELTA_ALPHA_I,
        const float &DELTA_ALPHA_J,
        const Vector< float> &Q_I,
        const Vector< float> &Q_J
      );

      void oclUpdateGradientBar
      (
        const float &EFFECTIVE_COST,
        const Vector< float> &Q_I
      );

      void oclUpdateGradientWithBias();

      void oclAddAlphaKernelToGradient
      (
        const Vector< float> &Q_I,
        const float &ALPHA
      );

      void oclAssembleFinalAlphaVector
      (
        Vector< float> &FINAL_ALPHAS,
        Vector< cl_uint> &SV_INDECES
      );

      void oclInitializeGradient
      (
        Vector< float> &Q_I
      );

      void oclInitializeGradientBar
      (
        Vector< float> &Q_I
      );

      //! @brief responsible for updating to a valid queue
      //! @param TOOLS opencl tools
      void UpdateQueue( Tools &TOOLS);

      //! @brief Set the members of this object from the given LABEL
      //! @param LABEL the label containing members that should be read of this class
      //! @param ERROR_STREAM stream with which to write errors
      //! @return result of any validation performed internally
      io::ValidationResult PreReadHook( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM)
      {
        UpdateQueue( GetTools());
        return io::ValidationResult( true);
      }

    }; // class ApproximatorSequentialMinimialOptimization

  } // namespace opencl
} // namespace bcl

#endif // BCL_OPENCL_APPROXIMATOR_SEQUENTIAL_MINIMIAL_OPTIMIZATION_H_
