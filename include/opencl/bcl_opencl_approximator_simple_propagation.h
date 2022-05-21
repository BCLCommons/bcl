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

#ifndef BCL_OPENCL_APPROXIMATOR_SIMPLE_PROPAGATION_H_
#define BCL_OPENCL_APPROXIMATOR_SIMPLE_PROPAGATION_H_

// include the namespace header
#include "bcl_opencl.h"

// include other forward headers - sorted alphabetically
#include "math/bcl_math.fwd.hh"
#include "model/bcl_model.fwd.hh"
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_opencl_command_queue.h"
#include "bcl_opencl_matrix.h"
#include "bcl_opencl_matrix_add.h"
#include "bcl_opencl_matrix_multiply.h"
#include "bcl_opencl_matrix_transpose.h"
#include "bcl_opencl_model_interface.h"
#include "bcl_opencl_rmsd.h"
#include "bcl_opencl_transfer_function_sigmoid.h"
#include "bcl_opencl_vector.h"
#include "bcl_opencl_vector_matrix_add.h"
#include "descriptor/bcl_descriptor_dataset.h"
#include "model/bcl_model_approximator_base.h"
#include "model/bcl_model_feature_data_set.h"
#include "model/bcl_model_objective_function_wrapper.h"
#include "model/bcl_model_transfer_function_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opencl
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ApproximatorSimplePropagation
    //! @brief opencl implementation of simple propagation for training neural networks - optimized for GPU
    //!
    //! @see @link example_opencl_approximator_simple_propagation.cpp @endlink
    //! @author loweew
    //! @date Mar 28, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ApproximatorSimplePropagation :
      public model::ApproximatorBase
    {
    private:

    //////////
    // data //
    //////////

      //! hidden architecture
      storage::Vector< size_t>                            m_HiddenArchitecture;

      //! alpha
      float                                               m_Alpha;

      //! eta
      float                                               m_Eta;

      //! steps per call
      size_t                                              m_StepsPerCall;

      //! training features buffer
      Matrix< float>                                      m_TrainingFeaturesOnDevice;

      //! training target buffer
      Matrix< float>                                      m_TrainingTargetOnDevice;

      //! monitor features buffer
      Matrix< float>                                      m_MonitorFeaturesOnDevice;

      //! monitor target buffer
      Matrix< float>                                      m_MonitorTargetOnDevice;

      //! the weight matrixes
      mutable storage::Vector< linal::Matrix< float> >    m_Weight;

      //! the bias
      mutable storage::Vector< linal::Vector< float> >    m_Bias;

      //! the weight matrixes
      storage::Vector< Matrix< float> >                   m_WeightBuffers;

      //! the transposed weight matrixes
      storage::Vector< Matrix< float> >                   m_TransWeightBuffers;

      //! the bias
      storage::Vector< Vector< float> >                   m_BiasBuffers;

      //! hidden buffers
      storage::Vector< Matrix< float> >                   m_HiddenBuffers;

      //! hidden input buffers
      storage::Vector< Matrix< float> >                   m_HiddenInputBuffers;

      //! errors buffers
      storage::Vector< Matrix< float> >                   m_ErrorsBuffers;

      //! transposed errors buffers
      storage::Vector< Matrix< float> >                   m_TransErrorsBuffers;

      //! slopes bias buffers
      storage::Vector< Vector< float> >                   m_SlopesBiasBuffers;

      //! change bias buffers
      storage::Vector< Vector< float> >                   m_ChangeBiasBuffers;

      //! slopes weight buffers
      storage::Vector< Matrix< float> >                   m_SlopesWeightBuffers;

      //! change weight buffers
      storage::Vector< Matrix< float> >                   m_ChangeWeightBuffers;

      //! opencl queue
      CommandQueue                                        m_Queue;

      //! gpu rmsd calc
      RMSD                                                m_GpuRMSD;

      //! gpu mmult
      MatrixMultiply< float>                              m_GpuMMult;

      //! gpu transpose
      MatrixTranspose< float>                             m_GpuTranspose;

      //! gpu m-m add
      MatrixAdd< float>                                   m_GpuMMAdd;

      //! gpu v-m add
      VectorMatrixAdd< float>                             m_GpuVMAdd;

      //! gpu sigmoid
      TransferFunctionSigmoid< float>                     m_GpuSigmoid;

      //! opencl program
      cl::Program                                         m_Program;

      //! network filename
      std::string                                         m_NetworkFilename;

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

      //! block size for kernels
      static const cl_uint s_Blocksize;

    public:

    //////////
    // data //
    //////////

      //! opencl compiler flags
      static const char *s_CLCompilerOptions;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      ApproximatorSimplePropagation();

      //! @brief constructor from training data, transfer, rescale, and objective functions
      //! @param TRAINING_DATA data to train the NeuralNetwork on
      //! @param ARCHITECTURE the architecture of the neural network, determines size of the matrices
      //! @param ALPHA learning momentum
      //! @param ETA learning rate
      //! @param STEPS_PER_CALL steps per call
      //! @param QUEUE command queue
      ApproximatorSimplePropagation
      (
        util::ShPtr< descriptor::Dataset> &TRAINING_DATA,
        const storage::Vector< size_t> &ARCHITECTURE,
        const float ALPHA,
        const float ETA,
        const size_t STEPS_PER_CALL,
        const CommandQueue &QUEUE
      );

      //! @brief copy constructor
      //! @return a new ApproximatorSimplePropagation copied from this instance
      ApproximatorSimplePropagation *Clone() const;

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
      void SetObjectiveFunction
      (
        const util::ShPtr< model::ObjectiveFunctionWrapper> &OBJ
      )
      {
        model::ApproximatorBase::SetObjectiveFunction( OBJ);
        const size_t f_rows( m_ObjectiveFunction->GetData()->GetFeaturesPtr()->GetNumberFeatures());
        const size_t f_cols( m_ObjectiveFunction->GetData()->GetFeaturesPtr()->GetFeatureSize());
        const size_t f_row_pad( ( s_Blocksize - ( f_rows % s_Blocksize)) % s_Blocksize);
        const size_t f_col_pad( ( s_Blocksize - ( f_cols % s_Blocksize)) % s_Blocksize);

        const size_t r_rows( m_ObjectiveFunction->GetData()->GetResultsPtr()->GetNumberFeatures());
        const size_t r_cols( m_ObjectiveFunction->GetData()->GetResultsPtr()->GetFeatureSize());
        const size_t r_row_pad( ( s_Blocksize - ( r_rows % s_Blocksize)) % s_Blocksize);
        const size_t r_col_pad( ( s_Blocksize - ( r_cols % s_Blocksize)) % s_Blocksize);

        m_MonitorFeaturesOnDevice = Matrix< float>( m_ObjectiveFunction->GetData()->GetFeaturesPtr()->GetMatrix(), m_Queue, f_row_pad, f_col_pad);
        m_MonitorTargetOnDevice = Matrix< float>( m_ObjectiveFunction->GetData()->GetResultsPtr()->GetMatrix(), m_Queue, r_row_pad, r_col_pad);
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief construct a model from the current iterate
      //! @return shptr to the new model interface
      util::ShPtr< model::Interface> GetCurrentModel() const;

      //! @brief construct a model from the current iterate
      //! @return shptr to the new model interface
      util::ShPtr< model::Interface> GetCurrentGPUModel() const;

      //! @brief returns the current approximation
      //! @return current argument result pair
      const util::ShPtr< storage::Pair< util::ShPtr< model::Interface>, float> > GetCurrentApproximation() const;

      //! @brief conducts the next approximation step and stores the approximation
      void Next();

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read NeuralNetwork from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write NeuralNetwork into std::ostream
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
      bool CanContinue() const
      {
        return true;
      }

      //! run forward through ANN and compute err terms m_Hidden (test)
      void CalcHiddenTerms();

      //! run backward through ANN and compute err terms m_Errors (train)
      void CalcErrorTerms();

      //! compute changes m_Changes to be applied on ANN (train)
      void CalcChangeTerms();

      //! train ANN with a feature
      void Train();

      //! @brief update the weights and bias
      void UpdateWeights();

      //! @brief calculates the output layer error terms
      void OutputLayerErrorTerms();

      //! @brief calculates the hidden layer error terms
      //! @param LAYER the layer number
      void OtherLayerErrorTerms( const size_t &LAYER);

      //! @brief adds bias change terms to previous bias
      //! @param LAYER layer number
      void AddBiasChangeTerms( const size_t &LAYER);

      //! @brief helper function to do a column-wise reduction on the bias errors
      //! @param ERRORS the errors buffer
      //! @param REDUCED_ERRORS the reduced error vector
      void ReduceErrorsToVector( Matrix< float> &ERRORS, Vector< float> &REDUCED_ERRORS);

      //! @brief applies change terms
      //! @param LAYER layer number
      void ApplyWeightChanges( const size_t &LAYER);

      //! @brief applies change terms
      //! @param LAYER layer number
      void ApplyBiasChanges( const size_t &LAYER);

      //! @brief sets up the neural network with a particular architecture, after training data was set
      void SetupArchitecture();

      //! @brief Set the members of this object from the given LABEL
      //! @param LABEL the label containing members that should be read of this class
      //! @param ERROR_STREAM stream with which to write errors
      //! @return result of any validation performed internally
      io::ValidationResult PreReadHook( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM)
      {
        UpdateQueue( GetTools());
        return io::ValidationResult( true);
      }

      //! @brief responsible for updating to a valid queue
      //! @param TOOLS opencl tools
      void UpdateQueue( Tools &TOOLS);

    }; // class IterateSimpleProagation

  } // namespace opencl
} // namespace bcl

#endif // BCL_OPENCL_APPROXIMATOR_SIMPLE_PROPAGATION_H_
