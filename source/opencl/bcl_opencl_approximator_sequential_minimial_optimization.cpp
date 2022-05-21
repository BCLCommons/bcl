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

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "opencl/bcl_opencl_approximator_sequential_minimial_optimization.h"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_vector_const_interface.hpp"
#include "model/bcl_model_support_vector_kernel_rbf.h"
#include "model/bcl_model_support_vector_machine.h"
#include "opencl/bcl_opencl_kernel_sources.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opencl
  {
    const cl_uint ApproximatorSequentialMinimialOptimization::s_Lower_Bound( 0);    //!<  0 indicates lower boundary
    const cl_uint ApproximatorSequentialMinimialOptimization::s_Upper_Bound( 1);    //!<  1 indicates upper boundary
    const cl_uint ApproximatorSequentialMinimialOptimization::s_Free( 2);           //!<  2 indicates no boundary
    const float ApproximatorSequentialMinimialOptimization::m_EPS_A( 1e-12);
    const float ApproximatorSequentialMinimialOptimization::m_P( 0.1);
    const char *ApproximatorSequentialMinimialOptimization::s_CLCompilerOptions =
            "-cl-mad-enable -cl-fast-relaxed-math";

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> ApproximatorSequentialMinimialOptimization::s_Instance
    (
      util::Enumerated< model::ApproximatorBase>::AddInstance( new ApproximatorSequentialMinimialOptimization())
    );

    //! @brief default constructor
    ApproximatorSequentialMinimialOptimization::ApproximatorSequentialMinimialOptimization() :
      m_CostParameterC( 0.0),
      m_Status(),
      m_Alpha(),
      m_Gradient(),
      m_GradientBar(),
      m_Bias(),
      m_Signs(),
      m_Labels(),
      m_ActiveSize(),
      m_ProbLength(),
      m_OptimizationGapThreshold( 0.1),
      m_OptimizationGap( 1.0),
      m_Model( new SupportVectorMachine()),
      m_NumberIterations( 0),
      m_NumberCurrentSupportVectors( 0)
    {
    }

    //! @brief Iterate for Sequential Minimal Optimization Learning Algorithm
    //! @param COST_PARAMETER_C penalty parameter c for svm regression training
    //! @param MODEL initial support vector model
    //! @param TRAINING_DATA training data set of choice
    //! @param NUMBER_ITERATIONS
    ApproximatorSequentialMinimialOptimization::ApproximatorSequentialMinimialOptimization
    (
      const float COST_PARAMETER_C,
      const util::ShPtr< SupportVectorMachine> &MODEL,
      util::ShPtr< descriptor::Dataset> &TRAINING_DATA,
      const size_t NUMBER_ITERATIONS,
      const CommandQueue &QUEUE
    ) :
      m_CostParameterC( COST_PARAMETER_C),
      m_Status(),
      m_Alpha(),
      m_Gradient(),
      m_GradientBar(),
      m_Bias(),
      m_Signs(),
      m_Labels(),
      m_ActiveSize( 0),
      m_ProbLength( 2 * TRAINING_DATA->GetSize()),
      m_Queue( QUEUE),
      m_OptimizationGapThreshold( 0.1),
      m_OptimizationGap( 1.0),
      m_Model( MODEL),
      m_NumberIterations( NUMBER_ITERATIONS),
      m_NumberCurrentSupportVectors( 0)
    {
      cl_int error_number( CL_SUCCESS);
      m_Program = KernelSources::Compile( GetKernelSources().e_SequentialMinimalOptimization, util::CPPDataTypes::e_Float, m_Queue, std::string( s_CLCompilerOptions), &error_number);
      BCL_Assert( error_number == CL_SUCCESS, "error compiling programs:\n" + Tools::ErrorString( error_number));

      SetTrainingData( TRAINING_DATA);
    }

    //! @brief copy constructor
    //! @return a new ApproximatorSequentialMinimialOptimization copied from this instance
    ApproximatorSequentialMinimialOptimization *ApproximatorSequentialMinimialOptimization::Clone() const
    {
      return new ApproximatorSequentialMinimialOptimization( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ApproximatorSequentialMinimialOptimization::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &ApproximatorSequentialMinimialOptimization::GetAlias() const
    {
      static const std::string s_Name( "OpenclApproximatorSequentialMinimialOptimization");
      return s_Name;
    }

    //! @brief set training data set for a specific iterate in approximater framework
    //! @param DATA training data set
    void ApproximatorSequentialMinimialOptimization::SetTrainingData( util::ShPtr< descriptor::Dataset> &DATA)
    {
      m_TrainingData = DATA;
      DATA->GetFeatures().Rescale( model::SupportVectorMachine::s_DefaultInputRange);
      DATA->GetResults().Rescale( model::SupportVectorMachine::s_DefaultInputRange);

      // set rescale functions and kernel svm model used in iterative training process
      if( m_Model->GetNumberSupportVectors() == 0)
      {
        m_Model = util::ShPtr< SupportVectorMachine>
        (
          new SupportVectorMachine
          (
            0.0,
            linal::Vector< float>( 1),
            model::FeatureDataSet< float>( linal::Matrix< float>( 1, m_TrainingData->GetFeatureSize())),
            util::Implementation< model::SupportVectorKernelBase>( model::SupportVectorKernelRBF( m_Gamma)),
            *GetRescaleFeatureDataSet(),
            *GetRescaleResultDataSet(),
            m_Queue
          )
        );
      }

      m_TrainingFeaturesOnDevice = Matrix< float>( m_TrainingData->GetFeaturesPtr()->GetMatrix(), m_Queue);
      m_TrainingResultsOnDevice = Matrix< float>( m_TrainingData->GetResultsPtr()->GetMatrix(), m_Queue);

      // initialize model
      InitializeMemberVectorsForTraining();
    }

    //! @brief construct a model from the current iterate
    //! @return shptr to the new model interface
    util::ShPtr< model::Interface> ApproximatorSequentialMinimialOptimization::GetCurrentModel() const
    {
      return m_Model;
    }

    //! @brief returns the current approximation
    //! @return current argument result pair
    const util::ShPtr< storage::Pair< util::ShPtr< model::Interface>, float> >
      ApproximatorSequentialMinimialOptimization::GetCurrentApproximation() const
    {
      util::ShPtr< model::Interface> model( GetCurrentModel());
      return
        util::ShPtr< storage::Pair< util::ShPtr< model::Interface>, float> >
        (
          new storage::Pair< util::ShPtr< model::Interface>, float>
          (
            model.HardCopy(),
            m_OptimizationGap
          )
        );
    }

    //! @brief conducts the next approximation step and stores the approximation
    void ApproximatorSequentialMinimialOptimization::Next()
    {
      // optimization gap for indicating the error difference delta epsilon to a given epsilon
      // iterate for a number of internal iterations
      if( m_OptimizationGap > m_OptimizationGapThreshold && m_NumberCurrentSupportVectors < float( 0.9) * m_TrainingData->GetSize())
      {
        for( size_t counter( 0); counter < m_NumberIterations; ++counter)
        {
          m_OptimizationGap = IterationStep();
        }
      }

      // postprocess model and determine final support vectors
      FinalizeSupportVectorModel();

      m_NumberCurrentSupportVectors = m_Model->GetNumberSupportVectors();
      // create final pair with model and objective function evaluation
      util::ShPtr< storage::Pair< util::ShPtr< model::Interface>, float> > current_model
      (
        new storage::Pair< util::ShPtr< model::Interface>, float>( m_Model.HardCopy(), m_OptimizationGap)
      );
      this->GetTracker().Track( current_model);

      BCL_MessageStd
      (
        " #SV: " + util::Format()( m_Model->GetNumberSupportVectors())
        + " gap: " + util::Format()( m_OptimizationGap)
      );
    }

    //! @brief iterates one cycle and returns ShPtr to pair of resultant argument and corresponding score
    float ApproximatorSequentialMinimialOptimization::IterationStep()
    {
      size_t first_vector_index( 0);
      size_t second_vector_index( 0);

      // determine index i and j for the two data vectors for which the quadratic problem has to be solved
      const float optimization_gap
      (
        DetermineFeatureVectorCombination( first_vector_index, second_vector_index)
      );

      // Solving Quadratic Problem sub problems for two given data vectors
      SolveQuadraticProblemSubProblem( first_vector_index, second_vector_index);

      // return optimization difference to error threshold epsilon
      return optimization_gap;
    }

    //! @brief evaluates whether the approximation can continue
    //! @return true, if the approximation can continue - otherwise false
    bool ApproximatorSequentialMinimialOptimization::CanContinue() const
    {
      return m_OptimizationGap > m_OptimizationGapThreshold
             && m_NumberCurrentSupportVectors < float( 0.9) * m_TrainingData->GetSize();
    }

    //! @brief finalize support vector model and determine final SVs, Alphas and Bias
    void ApproximatorSequentialMinimialOptimization::FinalizeSupportVectorModel()
    {
      // reconstruct the whole gradient
      ReconstructGradient();

      // counter for iteration
      size_t vector_counter( 0);

      // final vector with all alphas of support vectors
      storage::Vector< float> alpha_final;

      storage::List< size_t> sv_indices;

      Vector< float> device_final_alphas( m_AlphaOnDevice.GetSize() / 2, m_Queue);
      Vector< cl_uint> device_sv_indices( m_AlphaOnDevice.GetSize() / 2, m_Queue);

      oclAssembleFinalAlphaVector( device_final_alphas, device_sv_indices);

      linal::Vector< float> tmp_alphas( device_final_alphas.GetHostVector());
      linal::Vector< cl_uint> tmp_sv_indices( device_sv_indices.GetHostVector());

      for( size_t count( 0), count_end( m_AlphaOnDevice.GetSize() / 2); count < count_end; ++count)
      {
        if( tmp_alphas( count) != 0)
        {
          alpha_final.PushBack( tmp_alphas( count));
          sv_indices.PushBack( tmp_sv_indices( count));
        }
      }

      // matrix containing support vectors
      linal::Matrix< float> sv_matrix
      (
        sv_indices.GetSize(),                             // rows
        m_TrainingData->GetFeatureSize(), // cols
        float( 0)                                  // default value
      );

      vector_counter = 0;

      // fill support vector matrix with vectors by support vector index
      for
      (
        storage::List< size_t>::const_iterator itr_sv( sv_indices.Begin()), itr_sv_end( sv_indices.End());
        itr_sv != itr_sv_end;
        ++itr_sv, ++vector_counter
      )
      {
        const model::FeatureReference< float> &support_vector( m_TrainingData->GetFeaturesPtr()->operator ()( *itr_sv));

        std::copy( support_vector.Begin(), support_vector.End(), sv_matrix[ vector_counter]);
      }

      // assemble SV Model
      m_Model->SetAlpha( linal::Vector< float>( alpha_final.Begin(), alpha_final.End()));
      m_Model->SetBias( CalculateBias());
      m_Model->SetSupportVectors( model::FeatureDataSet< float>( sv_matrix));
      m_Model->SetNumberSupportVectors( sv_matrix.GetNumberRows());
    }

    //! read NeuralNetwork from std::istream
    std::istream &ApproximatorSequentialMinimialOptimization::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_NumberIterations, ISTREAM);
      io::Serialize::Read( m_CostParameterC, ISTREAM);
      io::Serialize::Read( m_Status, ISTREAM);
      io::Serialize::Read( m_Alpha, ISTREAM);
      io::Serialize::Read( m_Gradient, ISTREAM);
      io::Serialize::Read( m_GradientBar, ISTREAM);
      io::Serialize::Read( m_Bias, ISTREAM);
      io::Serialize::Read( m_Signs, ISTREAM);
      io::Serialize::Read( m_Labels, ISTREAM);
      io::Serialize::Read( m_ActiveSize, ISTREAM);
      io::Serialize::Read( m_ProbLength, ISTREAM);
      io::Serialize::Read( m_Model, ISTREAM);
      io::Serialize::Read( m_TrainingData, ISTREAM);
      // return
      return ISTREAM;
    }

    //! write NeuralNetwork into std::ostream
    std::ostream &ApproximatorSequentialMinimialOptimization::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_NumberIterations, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_CostParameterC, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Status, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Alpha, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Gradient, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_GradientBar, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Bias, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Signs, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Labels, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_ActiveSize, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_ProbLength, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Model, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_TrainingData, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

    //! @brief method checks whether a LaGrange Multiplier reached a certain boundary or not
    //! @param ALPHA a LaGrange Multiplier
    //! @return const int that indicates whether ALPHA reached a certain boundary or not
    int ApproximatorSequentialMinimialOptimization::AlphaToStatus( const float &ALPHA) const
    {
      if( ALPHA >= m_CostParameterC - m_EPS_A)
      {
        return s_Upper_Bound;
      }
      else if( ALPHA <= m_EPS_A)
      {
        return s_Lower_Bound;
      }
      else
      {
        return s_Free;
      }
    }

    //! @brief
    void ApproximatorSequentialMinimialOptimization::UpdateAlphaStatus( const int &FEATURE_VECTOR_I)
    {
      if( m_Alpha( FEATURE_VECTOR_I) >= m_CostParameterC)
      {
        m_Status( FEATURE_VECTOR_I) = s_Upper_Bound;
      }
      else if( m_Alpha( FEATURE_VECTOR_I) <= 0)
      {
        m_Status( FEATURE_VECTOR_I) = s_Lower_Bound;
      }
      else
      {
        m_Status( FEATURE_VECTOR_I) = s_Free;
      }
    }

    //! @brief method checks whether a feature_vector i reached the upper boundary
    //! @param FEATURE_VECTOR_I a position of a feature vector in m_Status
    //! @return const bool that indicates whether FEATURE_VECTOR_I reached upper boundary
    inline bool ApproximatorSequentialMinimialOptimization::IsUpperBound( const size_t &FEATURE_VECTOR_I) const
    {
      return m_Status( FEATURE_VECTOR_I) == s_Upper_Bound;
    }

    //! @brief checks whether a feature_vector i reached the lower boundary
    //! @param FEATURE_VECTOR_I a position of a feature vector in m_Status
    //! @return const bool that indicates whether FEATURE_VECTOR_I reached lower boundary
    inline bool ApproximatorSequentialMinimialOptimization::IsLowerBound( const size_t &FEATURE_VECTOR_I) const
    {
      return m_Status( FEATURE_VECTOR_I) == s_Lower_Bound;
    }

    //! @brief checks whether a feature_vector i is not bound and between the boundaries
    //! @param FEATURE_VECTOR_I a position of a feature vector in m_Status
    //! @return const bool that indicates whether FEATURE_VECTOR_I reached no boundary
    inline bool ApproximatorSequentialMinimialOptimization::IsFree( const size_t &FEATURE_VECTOR_I) const
    {
      return m_Status( FEATURE_VECTOR_I) == s_Free;
    }

    //! @brief Initialize member vectors and variables to prepare SVR training
    //! @param TRAINING_DATA feature vector data set with labels
    //! @param SVR_MODEL model to be trained
    void ApproximatorSequentialMinimialOptimization::InitializeMemberVectorsForTraining()
    {
      // Initialization of variables and vectors
      const size_t number_feature_vectors( m_TrainingData->GetSize());

      // has to be set in case of this class is read from file and constructor was not explicitly applied
      m_ActiveSize = 0;
      m_ProbLength = 2 * number_feature_vectors;

      // initialize vector for lagrange multipliers
      m_Alpha = linal::Vector< float>( m_ProbLength, 0.0);

      // initialize vector for labeling of
      m_Labels = linal::Vector< float>( m_ProbLength, float( 1.0));
      m_Signs = linal::Vector< cl_int>( m_ProbLength, 1);

      for( size_t count( number_feature_vectors), count_end( m_ProbLength); count < count_end; ++count)
      {
        m_Signs( count) = -1;
        m_Labels( count) = -1.0;
      }

      // status indication of KKT determination for every vector
      m_Status = linal::Vector< cl_uint>( m_ProbLength);
      m_GradientBar = linal::Vector< float>( m_ProbLength, 0.0);

      // initialization of bias for support vector model
      m_Bias = linal::Vector< float>( m_ProbLength, m_P);

      // initializing
      for( size_t progress( 0); progress < number_feature_vectors; ++progress)
      {
        const float result( m_TrainingData->GetResultsPtr()->operator()( progress)( 0));
        m_Bias( progress) -= result;
        m_Bias( number_feature_vectors + progress) += result;
      }

      m_Gradient = m_Bias;

      m_Q_i                 = Vector< float>( m_ProbLength, m_Queue);
      m_Q_d                 = Vector< float>( m_ProbLength, m_Queue, 0, 1.0);
      m_Q_j                 = Vector< float>( m_ProbLength, m_Queue);
      m_GradientOnDevice    = Vector< float>( m_Gradient, m_Queue);
      m_BiasOnDevice        = Vector< float>( m_Bias, m_Queue);
      m_AlphaOnDevice       = Vector< float>( m_Alpha, m_Queue);
      m_LabelsOnDevice      = Vector< float>( m_Labels, m_Queue);
      m_SignsOnDevice       = Vector< cl_int>( m_Signs, m_Queue);
      m_GradientBarOnDevice = Vector< float>( m_GradientBar, m_Queue);
      m_StatusOnDevice      = Vector< cl_uint>( m_Status, m_Queue);

      // update alpha_status
      oclUpdateAllAlphaStatus();

      m_ActiveSize = m_ProbLength;

      Vector< float> kernel_vector( m_ProbLength / 2, m_Queue);

      // initializing Gradients
      for( size_t progress( 0); progress < m_ActiveSize; ++progress)
      {
        if( !IsLowerBound( progress))
        {
          oclGetInputIKernelMatrixResultingVector( progress, m_Gamma, m_Q_i, kernel_vector);

          oclInitializeGradient( m_Q_i);

          if( IsUpperBound( progress))
          {
            oclInitializeGradientBar( m_Q_i);
          }
        }
      }

      BCL_MessageDbg( " Initialization.. complete");
    }

    //! @brief determine the bias value for Support Vector Regression Model
    //! @return calculated bias value
    float ApproximatorSequentialMinimialOptimization::CalculateBias()
    {
      float bias;

      const cl_uint block_size( 128);
      const cl_uint num_groups( ( m_Alpha.GetSize() % block_size == 0 ? 0 : 1) + ( m_Alpha.GetSize() / block_size));
      Vector< float> upper_output( num_groups, m_Queue);
      Vector< float> lower_output( num_groups, m_Queue);
      Vector< cl_uint> nr_free_output( num_groups, m_Queue);
      Vector< float> sum_free_output( num_groups, m_Queue);
      float final_upper, final_lower, final_sum_free;
      cl_uint final_nr_free;
      oclCalculateBias( upper_output, lower_output, nr_free_output, sum_free_output, final_upper, final_lower, final_nr_free, final_sum_free);

      // calculate bias
      if( final_nr_free > 0)
      {
        bias = final_sum_free / final_nr_free;
      }
      else
      {
        bias = ( final_upper + final_lower) / 2;
      }

      // return calculated bias value
      return bias;
    }

    //! @brief heuristic approach to find a feature_vector i and j
    //! @brief depending on m_Gradients for examination of SMO Classification
    //! @param TRAINING_DATA data set of feature vectors inclusive labels
    //! @param SVR_MODEL
    //! @param FIRST_VECTOR_INDEX
    //! @param SECOND_VECTOR_INDEX
    //! @return pair of two indices for feature_vector i and j
    float ApproximatorSequentialMinimialOptimization::DetermineFeatureVectorCombination
    (
      size_t &FIRST_VECTOR_INDEX,
      size_t &SECOND_VECTOR_INDEX
    )
    {
      // return i,j such that
      // i: maximizes -y_i * grad(f)_i, i in I_up(\alpha)
      // j: minimizes the decrease of obj value
      //    (if quadratic coefficient <= 0, replace it with tau)
      //    -y_j*grad(f)_j < -y_i*grad(f)_i, j in I_low(\alpha)

      // input for kernel as std in opencl doesn't work
      float maximum_gradient( -std::numeric_limits< float>::infinity());
      cl_int maximum_gradient_index( -1);
      cl_int gmin_ind( -1);
      float gmax2( -std::numeric_limits< float>::infinity());

      oclFindMaxGradient( maximum_gradient, maximum_gradient_index);

      const cl_int index_i( maximum_gradient_index);

      Vector< float> kernel_vector( m_ProbLength / 2, m_Queue);

      if( index_i != -1) // null Q_i not accessed: Gmax=-INF if index_i=-1
      {
        oclGetInputIKernelMatrixResultingVector( index_i, m_Gamma, m_Q_i, kernel_vector);
      }

      const cl_uint block_size( 128);
      const cl_uint num_groups( ( m_Labels.GetSize() % block_size == 0 ? 0 : 1) + ( m_Labels.GetSize() / block_size));
      Vector< float> tmp_gmax( num_groups, m_Queue);
      Vector< float> tmp_obj_diff_min( num_groups, m_Queue);
      Vector< cl_uint> tmp_gmin_ind( num_groups, m_Queue);

      oclGetGminGmax2( m_Q_i, m_Q_d, index_i, maximum_gradient, tmp_gmax, tmp_obj_diff_min, tmp_gmin_ind, gmax2, gmin_ind);

      FIRST_VECTOR_INDEX  = maximum_gradient_index;
      SECOND_VECTOR_INDEX = gmin_ind;

      return maximum_gradient + gmax2; // optimization gap
    }

    //! @brief solve sub problem of quadratic problem according two feature vectors
    //! @param SVR_MODEL support vector machine model of interest
    //! @param FEATURE_VECTOR_I index of feature vector in MATRIX
    //! @param FEATURE_VECTOR_J index of feature vector in MATRIX
    //! @return a pair of vector of float with computed values of kernel function
    void ApproximatorSequentialMinimialOptimization::SolveQuadraticProblemSubProblem
    (
      const int &FEATURE_VECTOR_I,
      const int &FEATURE_VECTOR_J
    )
    {
      Vector< float> kernel_vector_i( m_Alpha.GetSize() / 2, m_Queue);

      oclGetInputIKernelMatrixResultingVector( FEATURE_VECTOR_I, m_Gamma, m_Q_i, kernel_vector_i);

      Vector< float> kernel_vector_j( m_Alpha.GetSize() / 2, m_Queue);
      oclGetInputIKernelMatrixResultingVector( FEATURE_VECTOR_J, m_Gamma, m_Q_j, kernel_vector_j);

      float quad_coef( m_Q_i.GetHostVector()( FEATURE_VECTOR_I) + m_Q_j.GetHostVector()( FEATURE_VECTOR_J) + 2 * m_Q_i.GetHostVector()( FEATURE_VECTOR_J));

      if( quad_coef <= 0)
      {
        quad_coef = m_EPS_A;
      }

      m_Gradient = m_GradientOnDevice.GetHostVector();
      m_Alpha    = m_AlphaOnDevice.GetHostVector();
      m_Labels   = m_LabelsOnDevice.GetHostVector();

      const float old_m_alpha_i( m_Alpha( FEATURE_VECTOR_I));
      const float old_m_alpha_j( m_Alpha( FEATURE_VECTOR_J));

      // if one of the two classLabels is negative
      // computing LaGrange Multipliers m_Alpha
      if( m_Labels( FEATURE_VECTOR_I) * m_Labels( FEATURE_VECTOR_J) < 0)
      {
        const float delta( ( -1 * m_Gradient( FEATURE_VECTOR_I) - m_Gradient( FEATURE_VECTOR_J)) / quad_coef);

        const float diff( m_Alpha( FEATURE_VECTOR_I) - m_Alpha( FEATURE_VECTOR_J));

        m_Alpha( FEATURE_VECTOR_I) += delta;
        m_Alpha( FEATURE_VECTOR_J) += delta;

        if( diff > 0)
        {
          if( m_Alpha( FEATURE_VECTOR_J) < 0)
          {
            m_Alpha( FEATURE_VECTOR_J) = 0;
            m_Alpha( FEATURE_VECTOR_I) = diff;
          }
          if( m_Alpha( FEATURE_VECTOR_I) > m_CostParameterC)
          {
            m_Alpha( FEATURE_VECTOR_I) = m_CostParameterC;
            m_Alpha( FEATURE_VECTOR_J) = m_CostParameterC - diff;
          }
        }
        else
        {
          if( m_Alpha( FEATURE_VECTOR_I) < 0)
          {
            m_Alpha( FEATURE_VECTOR_I) = 0;
            m_Alpha( FEATURE_VECTOR_J) = -diff;
          }
          if( m_Alpha( FEATURE_VECTOR_J) > m_CostParameterC)
          {
            m_Alpha( FEATURE_VECTOR_J) = m_CostParameterC;
            m_Alpha( FEATURE_VECTOR_I) = m_CostParameterC + diff;
          }
        }
      }
      else // if both classLabels are positive
      {
        const float delta( ( m_Gradient( FEATURE_VECTOR_I) - m_Gradient( FEATURE_VECTOR_J)) / quad_coef);
        const float sum( m_Alpha( FEATURE_VECTOR_I) + m_Alpha( FEATURE_VECTOR_J));

        m_Alpha( FEATURE_VECTOR_I) -= delta;
        m_Alpha( FEATURE_VECTOR_J) += delta;

        if( sum > m_CostParameterC)
        {
          if( m_Alpha( FEATURE_VECTOR_I) > m_CostParameterC)
          {
            m_Alpha( FEATURE_VECTOR_I) = m_CostParameterC;
            m_Alpha( FEATURE_VECTOR_J) = sum - m_CostParameterC;
          }
          if( m_Alpha( FEATURE_VECTOR_J) > m_CostParameterC)
          {
            m_Alpha( FEATURE_VECTOR_J) = m_CostParameterC;
            m_Alpha( FEATURE_VECTOR_I) = sum - m_CostParameterC;
          }
        }
        else
        {
          if( m_Alpha( FEATURE_VECTOR_J) < 0)
          {
            m_Alpha( FEATURE_VECTOR_J) = 0;
            m_Alpha( FEATURE_VECTOR_I) = sum;
          }
          if( m_Alpha( FEATURE_VECTOR_I) < 0)
          {
            m_Alpha( FEATURE_VECTOR_I) = 0;
            m_Alpha( FEATURE_VECTOR_J) = sum;
          }
        }
      }

      // update Gradient
      const float delta_alpha_i( m_Alpha( FEATURE_VECTOR_I) - old_m_alpha_i);
      const float delta_alpha_j( m_Alpha( FEATURE_VECTOR_J) - old_m_alpha_j);

      oclUpdateGradient( delta_alpha_i, delta_alpha_j, m_Q_i, m_Q_j);

      // update alpha_status and m_GradientBar
      m_Status = m_StatusOnDevice.GetHostVector();
      const bool feature_i_is_upper_bound( IsUpperBound( FEATURE_VECTOR_I));
      const bool feature_j_is_upper_bound( IsUpperBound( FEATURE_VECTOR_J));

      UpdateAlphaStatus( FEATURE_VECTOR_I);
      UpdateAlphaStatus( FEATURE_VECTOR_J);

      m_StatusOnDevice = Vector< cl_uint>( m_Status, m_Queue);
      m_AlphaOnDevice = Vector< float>( m_Alpha, m_Queue);

      short cost_multiplier( 0);

      if( feature_i_is_upper_bound != IsUpperBound( FEATURE_VECTOR_I))
      {
        cost_multiplier += ( feature_i_is_upper_bound ? -1 : 1);
      }
      if( feature_j_is_upper_bound != IsUpperBound( FEATURE_VECTOR_J))
      {
        cost_multiplier += ( feature_j_is_upper_bound ? -1 : 1);
      }

      if( cost_multiplier != 0)
      {
        const double effective_cost( cost_multiplier * m_CostParameterC);

        oclUpdateGradientBar( effective_cost, m_Q_i);
      }
    }

    void ApproximatorSequentialMinimialOptimization::ReconstructGradient()
    {
      // reconstruct inactive elements of m_Gradient from m_GradientBar and free variables
      if( m_ActiveSize == m_ProbLength)
      {
        return;
      }

      oclUpdateGradientWithBias();

      // vector of kernel values of vector i and all other vectors in training set
      Vector< float> kernel_vector( m_Alpha.GetSize() / 2, m_Queue);

      // progress counter accessing elements at a certain position in vector
      size_t progress( 0);

      // iterate over all LaGrange multipliers
      m_Alpha = m_AlphaOnDevice.GetHostVector();

      for
      (
        linal::Vector< float>::iterator
          iter_begin_alpha( m_Alpha.Begin()),
          iter_end_alpha( m_Alpha.End());
        iter_begin_alpha != iter_end_alpha;
        ++iter_begin_alpha, ++progress
      )
      {
        // if correspondent status values indicates that it is not between upper or lower bound
        if( IsFree( progress))
        {
          oclGetInputIKernelMatrixResultingVector( progress, m_Gamma, m_Q_i, kernel_vector);

          oclAddAlphaKernelToGradient( m_Q_i, *iter_begin_alpha);
        }
      }
    }

    void ApproximatorSequentialMinimialOptimization::oclFindMaxGradient( float &MAX_GRAD, cl_int &MAX_IND)
    {
      // catch errors
      cl_int error_number = CL_SUCCESS;

      // dimensions
      const cl_uint elements( m_GradientOnDevice.GetSize());
      const cl_uint block_size( 128);
      const size_t num_groups( ( elements % block_size == 0 ? 0 : 1) + ( elements / block_size));

      Vector< float> tmp_max( num_groups, m_Queue);
      Vector< cl_uint> tmp_indexes( num_groups, m_Queue);

      // create kernel
      cl::Kernel kernel( m_Program, "FindMaxGradient", &error_number);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + opencl::Tools::ErrorString( error_number));

      // set thread block dimensions
      const cl::NDRange block_dimensions( block_size); //< all thread blocks have same dimensions
      const cl::NDRange offset;
      const cl::NDRange worksize( Tools::RoundUp( block_size, elements));

      // add change in bias to previous bias
      error_number  = kernel.setArg( 0, tmp_max.GetData());
      error_number |= kernel.setArg( 1, tmp_indexes.GetData());
      error_number |= kernel.setArg( 2, m_GradientOnDevice.GetData());
      error_number |= kernel.setArg( 3, m_LabelsOnDevice.GetData());
      error_number |= kernel.setArg( 4, m_StatusOnDevice.GetData());
      error_number |= kernel.setArg( 5, s_Upper_Bound);
      error_number |= kernel.setArg( 6, s_Lower_Bound);
      error_number |= kernel.setArg( 7, block_size * sizeof( float), 0);
      error_number |= kernel.setArg( 8, block_size * sizeof( cl_uint), 0);
      error_number |= kernel.setArg( 9, elements);
      BCL_Assert( error_number == CL_SUCCESS, "FindMaxGradient error: " + opencl::Tools::ErrorString( error_number));

      error_number = m_Queue.enqueueNDRangeKernel( kernel, offset, worksize, block_dimensions, NULL, NULL);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + opencl::Tools::ErrorString( error_number));

      linal::Vector< float> max_vector( tmp_max.GetHostVector());
      linal::Vector< cl_uint> ind_vector( tmp_indexes.GetHostVector());

      // complete on cpu
      float max_element( -std::numeric_limits< float>::infinity());
      size_t final_index( 0);
      for( size_t count( 0); count < num_groups; ++count)
      {
        size_t greater_than( max_vector( count) >= max_element ? 1 : 0);
        greater_than ? max_element = max_vector( count), final_index = ind_vector( count) : 0;
      }
      MAX_GRAD = max_element;
      MAX_IND  = final_index;
      m_Queue.finish();
    }

    void ApproximatorSequentialMinimialOptimization::oclGetInputIKernelMatrixResultingVector
    (
      const cl_uint &VECTOR_ID,
      const float &GAMMA,
      Vector< float> &OUTPUT,
      Vector< float> &KERNEL_VECTOR
    )
    {
      // catch errors
      cl_int error_number = CL_SUCCESS;

      // dimensions
      const cl_uint rows( m_TrainingFeaturesOnDevice.GetNumberRows());
      const cl_uint cols( m_TrainingFeaturesOnDevice.GetNumberCols());

      const cl_uint block_size( 128);

      // create kernel
      cl::Kernel kernel_a( m_Program, "ComputeInputIKernelMatrix", &error_number);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + opencl::Tools::ErrorString( error_number));
      cl::Kernel kernel_b( m_Program, "GetInputIKernelMatrixResultingVector", &error_number);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + opencl::Tools::ErrorString( error_number));

      cl_uint index_vector;

      // adjust index for input vector i
      if( VECTOR_ID >= rows)
      {
        index_vector = VECTOR_ID - ( m_ProbLength / 2);
      }
      else
      {
        index_vector = VECTOR_ID;
      }

      // set thread block dimensions
      const cl::NDRange block_dimensions( block_size); //< all thread blocks have same dimensions
      const cl::NDRange offset;
      const cl::NDRange worksize_a( Tools::RoundUp( block_size, rows));
      const cl::NDRange worksize_b( Tools::RoundUp( block_size, m_ProbLength));

      // add change in bias to previous bias
      error_number  = kernel_a.setArg( 0, index_vector);
      error_number |= kernel_a.setArg( 1, m_TrainingFeaturesOnDevice.GetData());
      error_number |= kernel_a.setArg( 2, cols);
      error_number |= kernel_a.setArg( 3, rows);
      error_number |= kernel_a.setArg( 4, KERNEL_VECTOR.GetData());
      error_number |= kernel_a.setArg( 5, block_size * sizeof( float), 0);
      error_number |= kernel_a.setArg( 6, GAMMA);
      BCL_Assert( error_number == CL_SUCCESS, "ComputeInputIKernelMatrix error: " + opencl::Tools::ErrorString( error_number));

      // add change in bias to previous bias
      error_number  = kernel_b.setArg( 0, VECTOR_ID);
      error_number |= kernel_b.setArg( 1, KERNEL_VECTOR.GetData());
      error_number |= kernel_b.setArg( 2, m_SignsOnDevice.GetData());
      error_number |= kernel_b.setArg( 3, cl_uint( m_ProbLength));
      error_number |= kernel_b.setArg( 4, OUTPUT.GetData());
      BCL_Assert( error_number == CL_SUCCESS, "GetInputIKernelMatrixResultingVector error: " + opencl::Tools::ErrorString( error_number));

      error_number = m_Queue.enqueueNDRangeKernel( kernel_a, offset, worksize_a, block_dimensions, NULL, NULL);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + opencl::Tools::ErrorString( error_number));

      error_number = m_Queue.enqueueNDRangeKernel( kernel_b, offset, worksize_b, block_dimensions, NULL, NULL);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + opencl::Tools::ErrorString( error_number));
      m_Queue.finish();
    }

    void ApproximatorSequentialMinimialOptimization::oclGetGminGmax2
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
    )
    {
      // catch errors
      cl_int error_number = CL_SUCCESS;

      // dimensions
      const cl_uint elements( m_GradientOnDevice.GetSize());
      const cl_uint block_size( 128);
      const size_t num_groups( ( elements % block_size == 0 ? 0 : 1) + ( elements / block_size));

      // create kernel
      cl::Kernel kernel( m_Program, "GetGminGmax2", &error_number);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + opencl::Tools::ErrorString( error_number));

      // set thread block dimensions
      const cl::NDRange block_dimensions( block_size); //< all thread blocks have same dimensions
      const cl::NDRange offset;
      const cl::NDRange worksize( Tools::RoundUp( block_size, elements));

      // add change in bias to previous bias
      error_number  = kernel.setArg( 0, Q_I.GetData());
      error_number |= kernel.setArg( 1, Q_D.GetData());
      error_number |= kernel.setArg( 2, m_GradientOnDevice.GetData());
      error_number |= kernel.setArg( 3, m_LabelsOnDevice.GetData());
      error_number |= kernel.setArg( 4, m_StatusOnDevice.GetData());
      error_number |= kernel.setArg( 5, s_Upper_Bound);
      error_number |= kernel.setArg( 6, s_Lower_Bound);
      error_number |= kernel.setArg( 7, cl_uint( INDEX_I));
      error_number |= kernel.setArg( 8, MAX_GRADIENT);
      error_number |= kernel.setArg( 9, m_EPS_A);
      error_number |= kernel.setArg( 10, elements);
      error_number |= kernel.setArg( 11, TMP_GMAX2_OUTPUT.GetData());
      error_number |= kernel.setArg( 12, TMP_OBJ_DIFF_MIN_OUTPUT.GetData());
      error_number |= kernel.setArg( 13, TMP_GMIN_INDEX_OUTPUT.GetData());
      error_number |= kernel.setArg( 14, block_size * sizeof( float), 0);
      error_number |= kernel.setArg( 15, block_size * sizeof( float), 0);
      error_number |= kernel.setArg( 16, block_size * sizeof( cl_uint), 0);
      BCL_Assert( error_number == CL_SUCCESS, "GetGminGmax2 error: " + opencl::Tools::ErrorString( error_number));

      error_number = m_Queue.enqueueNDRangeKernel( kernel, offset, worksize, block_dimensions, NULL, NULL);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + opencl::Tools::ErrorString( error_number));

      linal::Vector< float> gmax2_vector( TMP_GMAX2_OUTPUT.GetHostVector());
      linal::Vector< float> obj_diff_min_vector( TMP_OBJ_DIFF_MIN_OUTPUT.GetHostVector());
      linal::Vector< cl_uint> gmin_ind_vector( TMP_GMIN_INDEX_OUTPUT.GetHostVector());

      // complete on cpu
      float gmax_element( -std::numeric_limits< float>::infinity());
      float obj_diff_min( std::numeric_limits< float>::infinity());
      int gmin_ind( 0);
      for( size_t count( 0); count < num_groups; ++count)
      {
        if( gmax2_vector( count) > gmax_element)
        {
          gmax_element = gmax2_vector( count);
        }

        if( obj_diff_min_vector( count) < obj_diff_min)
        {
          obj_diff_min = obj_diff_min_vector( count);
          gmin_ind = gmin_ind_vector( count);
        }
      }

      FINAL_GMAX2 = gmax_element;
      FINAL_GMIN_IND = gmin_ind;

      m_Queue.finish();
    }

    void ApproximatorSequentialMinimialOptimization::oclCalculateBias
    (
      Vector< float> &UPPER_OUTPUT,
      Vector< float> &LOWER_OUTPUT,
      Vector< cl_uint>      &NR_FREE_OUTPUT,
      Vector< float> &SUM_FREE_OUTPUT,
      float &FINAL_UPPER_OUTPUT,
      float &FINAL_LOWER_OUTPUT,
      cl_uint      &FINAL_NR_FREE_OUTPUT,
      float &FINAL_SUM_FREE_OUTPUT
    )
    {
      // catch errors
      cl_int error_number = CL_SUCCESS;

      // dimensions
      const cl_uint elements( m_StatusOnDevice.GetSize());
      const cl_uint block_size( 128);
      const size_t num_groups( ( elements % block_size == 0 ? 0 : 1) + ( elements / block_size));

      // create kernel
      cl::Kernel kernel( m_Program, "CalculateBias", &error_number);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + opencl::Tools::ErrorString( error_number));

      // set thread block dimensions
      const cl::NDRange block_dimensions( block_size); //< all thread blocks have same dimensions
      const cl::NDRange offset;
      const cl::NDRange worksize( Tools::RoundUp( block_size, elements));

      // add change in bias to previous bias
      error_number  = kernel.setArg( 0, m_StatusOnDevice.GetData());
      error_number |= kernel.setArg( 1, m_GradientOnDevice.GetData());
      error_number |= kernel.setArg( 2, m_LabelsOnDevice.GetData());
      error_number |= kernel.setArg( 3, s_Lower_Bound);
      error_number |= kernel.setArg( 4, s_Upper_Bound);
      error_number |= kernel.setArg( 5, elements);
      error_number |= kernel.setArg( 6, block_size * sizeof( float), 0);
      error_number |= kernel.setArg( 7, block_size * sizeof( float), 0);
      error_number |= kernel.setArg( 8, block_size * sizeof( cl_uint), 0);
      error_number |= kernel.setArg( 9, block_size * sizeof( float), 0);
      error_number |= kernel.setArg( 10, UPPER_OUTPUT.GetData());
      error_number |= kernel.setArg( 11, LOWER_OUTPUT.GetData());
      error_number |= kernel.setArg( 12, NR_FREE_OUTPUT.GetData());
      error_number |= kernel.setArg( 13, SUM_FREE_OUTPUT.GetData());
      BCL_Assert( error_number == CL_SUCCESS, "CalculateBias error: " + opencl::Tools::ErrorString( error_number));

      error_number = m_Queue.enqueueNDRangeKernel( kernel, offset, worksize, block_dimensions, NULL, NULL);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + opencl::Tools::ErrorString( error_number));

      linal::Vector< cl_uint> nr_free_vector( NR_FREE_OUTPUT.GetHostVector());
      linal::Vector< float> sum_free_vector( SUM_FREE_OUTPUT.GetHostVector());
      linal::Vector< float> upper_vector( UPPER_OUTPUT.GetHostVector());
      linal::Vector< float> lower_vector( LOWER_OUTPUT.GetHostVector());

      FINAL_NR_FREE_OUTPUT = nr_free_vector.Sum();
      FINAL_SUM_FREE_OUTPUT = sum_free_vector.Sum();

      // complete on cpu
      float final_upper( std::numeric_limits< float>::infinity());
      float final_lower( -std::numeric_limits< float>::infinity());
      for( size_t count( 0); count < num_groups; ++count)
      {
        if( upper_vector( count) < final_upper)
        {
          final_upper = upper_vector( count);
        }

        if( lower_vector( count) > final_lower)
        {
          final_lower = lower_vector( count);
        }
      }

      FINAL_LOWER_OUTPUT = final_lower;
      FINAL_UPPER_OUTPUT = final_upper;

      m_Queue.finish();
    }

    void ApproximatorSequentialMinimialOptimization::oclUpdateAllAlphaStatus
    (
    )
    {
      // catch errors
      cl_int error_number = CL_SUCCESS;

      // dimensions
      const cl_uint elements( m_AlphaOnDevice.GetSize());

      const cl_uint block_size( 128);

      // create kernel
      cl::Kernel kernel( m_Program, "UpdateAllAlphaStatus", &error_number);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + opencl::Tools::ErrorString( error_number));

      // set thread block dimensions
      const cl::NDRange block_dimensions( block_size); //< all thread blocks have same dimensions
      const cl::NDRange offset;
      const cl::NDRange worksize( Tools::RoundUp( block_size, elements));

      // add change in bias to previous bias
      error_number |= kernel.setArg( 0, m_CostParameterC);
      error_number |= kernel.setArg( 1, s_Upper_Bound);
      error_number |= kernel.setArg( 2, s_Lower_Bound);
      error_number |= kernel.setArg( 3, s_Free);
      error_number |= kernel.setArg( 4, m_AlphaOnDevice.GetData());
      error_number |= kernel.setArg( 5, m_StatusOnDevice.GetData());
      error_number |= kernel.setArg( 6, elements);
      BCL_Assert( error_number == CL_SUCCESS, "UpdateAllAlphaStatus error: " + opencl::Tools::ErrorString( error_number));

      error_number = m_Queue.enqueueNDRangeKernel( kernel, offset, worksize, block_dimensions, NULL, NULL);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + opencl::Tools::ErrorString( error_number));
      m_Queue.finish();
    }

    void ApproximatorSequentialMinimialOptimization::oclUpdateGradient
    (
      const float &DELTA_ALPHA_I,
      const float &DELTA_ALPHA_J,
      const Vector< float> &Q_I,
      const Vector< float> &Q_J
    )
    {
      // catch errors
      cl_int error_number = CL_SUCCESS;

      // dimensions
      const cl_uint elements( m_GradientOnDevice.GetSize());

      const cl_uint block_size( 128);

      // create kernel
      cl::Kernel kernel( m_Program, "UpdateGradient", &error_number);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + opencl::Tools::ErrorString( error_number));

      // set thread block dimensions
      const cl::NDRange block_dimensions( block_size); //< all thread blocks have same dimensions
      const cl::NDRange offset;
      const cl::NDRange worksize( Tools::RoundUp( block_size, elements));

      // add change in bias to previous bias
      error_number  = kernel.setArg( 0, DELTA_ALPHA_I);
      error_number |= kernel.setArg( 1, DELTA_ALPHA_J);
      error_number |= kernel.setArg( 2, m_GradientOnDevice.GetData());
      error_number |= kernel.setArg( 3, Q_I.GetData());
      error_number |= kernel.setArg( 4, Q_J.GetData());
      error_number |= kernel.setArg( 5, elements);
      BCL_Assert( error_number == CL_SUCCESS, "UpdateGradient error: " + opencl::Tools::ErrorString( error_number));

      error_number = m_Queue.enqueueNDRangeKernel( kernel, offset, worksize, block_dimensions, NULL, NULL);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + opencl::Tools::ErrorString( error_number));
      m_Queue.finish();

    }

    void ApproximatorSequentialMinimialOptimization::oclUpdateGradientBar
    (
      const float &EFFECTIVE_COST,
      const Vector< float> &Q_I
    )
    {
      // catch errors
      cl_int error_number = CL_SUCCESS;

      // dimensions
      const cl_uint elements( m_GradientBarOnDevice.GetSize());

      const cl_uint block_size( 128);

      // create kernel
      cl::Kernel kernel( m_Program, "UpdateGradientBar", &error_number);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + opencl::Tools::ErrorString( error_number));

      // set thread block dimensions
      const cl::NDRange block_dimensions( block_size); //< all thread blocks have same dimensions
      const cl::NDRange offset;
      const cl::NDRange worksize( Tools::RoundUp( block_size, elements));

      // add change in bias to previous bias
      error_number  = kernel.setArg( 0, EFFECTIVE_COST);
      error_number |= kernel.setArg( 1, m_GradientBarOnDevice.GetData());
      error_number |= kernel.setArg( 2, Q_I.GetData());
      error_number |= kernel.setArg( 3, elements);
      BCL_Assert( error_number == CL_SUCCESS, "UpdateGradientBar error: " + opencl::Tools::ErrorString( error_number));

      error_number = m_Queue.enqueueNDRangeKernel( kernel, offset, worksize, block_dimensions, NULL, NULL);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + opencl::Tools::ErrorString( error_number));
      m_Queue.finish();
    }

    void ApproximatorSequentialMinimialOptimization::oclUpdateGradientWithBias()
    {
      // catch errors
      cl_int error_number = CL_SUCCESS;

      // dimensions
      const cl_uint elements( m_GradientBarOnDevice.GetSize());

      const cl_uint block_size( 128);

      // create kernel
      cl::Kernel kernel( m_Program, "UpdateGradientWithBias", &error_number);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + opencl::Tools::ErrorString( error_number));

      // set thread block dimensions
      const cl::NDRange block_dimensions( block_size); //< all thread blocks have same dimensions
      const cl::NDRange offset;
      const cl::NDRange worksize( Tools::RoundUp( block_size, elements));

      // add change in bias to previous bias
      error_number  = kernel.setArg( 0, m_GradientBarOnDevice.GetData());
      error_number |= kernel.setArg( 1, m_BiasOnDevice.GetData());
      error_number |= kernel.setArg( 2, elements);
      error_number |= kernel.setArg( 3, m_GradientOnDevice.GetData());
      BCL_Assert( error_number == CL_SUCCESS, "UpdateGradientBar error: " + opencl::Tools::ErrorString( error_number));

      error_number = m_Queue.enqueueNDRangeKernel( kernel, offset, worksize, block_dimensions, NULL, NULL);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + opencl::Tools::ErrorString( error_number));
      m_Queue.finish();
    }

    void ApproximatorSequentialMinimialOptimization::oclAddAlphaKernelToGradient
    (
      const Vector< float> &Q_I,
      const float &ALPHA
    )
    {
      // catch errors
      cl_int error_number = CL_SUCCESS;

      // dimensions
      const cl_uint elements( Q_I.GetSize());

      const cl_uint block_size( 128);

      // create kernel
      cl::Kernel kernel( m_Program, "AddAlphaKernelToGradient", &error_number);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + opencl::Tools::ErrorString( error_number));

      // set thread block dimensions
      const cl::NDRange block_dimensions( block_size); //< all thread blocks have same dimensions
      const cl::NDRange offset;
      const cl::NDRange worksize( Tools::RoundUp( block_size, elements));

      // add change in bias to previous bias
      error_number  = kernel.setArg( 0, cl_uint( m_ActiveSize));
      error_number |= kernel.setArg( 1, m_GradientOnDevice.GetData());
      error_number |= kernel.setArg( 2, Q_I.GetData());
      error_number |= kernel.setArg( 3, ALPHA);
      error_number |= kernel.setArg( 4, elements);
      BCL_Assert( error_number == CL_SUCCESS, "AddAlphaKernelToGradient error: " + opencl::Tools::ErrorString( error_number));

      error_number = m_Queue.enqueueNDRangeKernel( kernel, offset, worksize, block_dimensions, NULL, NULL);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + opencl::Tools::ErrorString( error_number));
      m_Queue.finish();
    }

    void ApproximatorSequentialMinimialOptimization::oclAssembleFinalAlphaVector
    (
      Vector< float> &FINAL_ALPHAS,
      Vector< cl_uint> &SV_INDECES
    )
    {
      // catch errors
      cl_int error_number = CL_SUCCESS;

      // dimensions
      const cl_uint elements( FINAL_ALPHAS.GetSize());

      const cl_uint block_size( 128);

      // create kernel
      cl::Kernel kernel( m_Program, "AssembleFinalAlphaVector", &error_number);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + opencl::Tools::ErrorString( error_number));

      // set thread block dimensions
      const cl::NDRange block_dimensions( block_size); //< all thread blocks have same dimensions
      const cl::NDRange offset;
      const cl::NDRange worksize( Tools::RoundUp( block_size, elements));

      // add change in bias to previous bias
      error_number  = kernel.setArg( 0, cl_uint( m_ProbLength / 2));
      error_number |= kernel.setArg( 1, m_AlphaOnDevice.GetData());
      error_number |= kernel.setArg( 2, FINAL_ALPHAS.GetData());
      error_number |= kernel.setArg( 3, SV_INDECES.GetData());
      BCL_Assert( error_number == CL_SUCCESS, "AssembleFinalAlphaVector error: " + opencl::Tools::ErrorString( error_number));

      error_number = m_Queue.enqueueNDRangeKernel( kernel, offset, worksize, block_dimensions, NULL, NULL);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + opencl::Tools::ErrorString( error_number));
      m_Queue.finish();
    }

    void ApproximatorSequentialMinimialOptimization::oclInitializeGradient
    (
      Vector< float> &Q_I
    )
    {
      // catch errors
      cl_int error_number = CL_SUCCESS;

      // dimensions
      const cl_uint elements( Q_I.GetSize());

      const cl_uint block_size( 128);

      // create kernel
      cl::Kernel kernel( m_Program, "InitializeGradient", &error_number);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + opencl::Tools::ErrorString( error_number));

      // set thread block dimensions
      const cl::NDRange block_dimensions( block_size); //< all thread blocks have same dimensions
      const cl::NDRange offset;
      const cl::NDRange worksize( Tools::RoundUp( block_size, elements));

      // add change in bias to previous bias
      error_number  = kernel.setArg( 0, m_AlphaOnDevice.GetData());
      error_number |= kernel.setArg( 1, Q_I.GetData());
      error_number |= kernel.setArg( 2, m_GradientOnDevice.GetData());
      error_number |= kernel.setArg( 3, elements);
      BCL_Assert( error_number == CL_SUCCESS, "InitializeGradient error: " + opencl::Tools::ErrorString( error_number));

      error_number = m_Queue.enqueueNDRangeKernel( kernel, offset, worksize, block_dimensions, NULL, NULL);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + opencl::Tools::ErrorString( error_number));
      m_Queue.finish();
    }

    void ApproximatorSequentialMinimialOptimization::oclInitializeGradientBar
    (
      Vector< float> &Q_I
    )
    {
      // catch errors
      cl_int error_number = CL_SUCCESS;

      // dimensions
      const cl_uint elements( Q_I.GetSize());

      const cl_uint block_size( 128);

      // create kernel
      cl::Kernel kernel( m_Program, "InitializeGradientBar", &error_number);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + opencl::Tools::ErrorString( error_number));

      // set thread block dimensions
      const cl::NDRange block_dimensions( block_size); //< all thread blocks have same dimensions
      const cl::NDRange offset;
      const cl::NDRange worksize( Tools::RoundUp( block_size, elements));

      // add change in bias to previous bias
      error_number  = kernel.setArg( 0, m_CostParameterC);
      error_number |= kernel.setArg( 1, Q_I.GetData());
      error_number |= kernel.setArg( 2, m_GradientBarOnDevice.GetData());
      error_number |= kernel.setArg( 3, elements);
      BCL_Assert( error_number == CL_SUCCESS, "InitializeGradientBar error: " + opencl::Tools::ErrorString( error_number));

      error_number = m_Queue.enqueueNDRangeKernel( kernel, offset, worksize, block_dimensions, NULL, NULL);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + opencl::Tools::ErrorString( error_number));
      m_Queue.finish();
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer ApproximatorSequentialMinimialOptimization::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "trains a support vector machine using sequential-minimal-optimization "
        "(see http://en.wikipedia.org/wiki/Sequential_Minimal_Optimization)"
      );

      parameters.AddInitializer
      (
        "objective function",
        "function that evaluates the model after each batch step",
        io::Serialization::GetAgent( &m_ObjectiveFunction->GetImplementation()),
        "RMSD"
      );
      parameters.AddInitializer
      (
        "gamma",
        "currently hardcoded to RBF, so this is the gamma value",
        io::Serialization::GetAgent( &m_Gamma),
        "0.5"
      );

      parameters.AddInitializer
      (
        "iterations",
        "# of iterations used internally to improve the optimization gap",
        io::Serialization::GetAgent( &m_NumberIterations),
        "0"
      );
      parameters.AddInitializer
      (
        "cost",
        "controls the trade off between allowing training errors and forcing rigid margins; high values may lead to better training values, but risks overtraining",
        io::Serialization::GetAgent( &m_CostParameterC),
        "0.0"
      );
      return parameters;
    }

    //! @brief responsible for updating to a valid queue
    //! @param TOOLS opencl tools
    void ApproximatorSequentialMinimialOptimization::UpdateQueue( Tools &TOOLS)
    {
      if( !TOOLS.HasCommandQueues())
      {
        return;
      }

      m_Queue = GetTools().GetFirstCommandQueue();

      cl_int error_number( CL_SUCCESS);
      m_Program = KernelSources::Compile( GetKernelSources().e_SequentialMinimalOptimization, util::CPPDataTypes::e_Float, m_Queue, std::string( s_CLCompilerOptions), &error_number);
      BCL_Assert( error_number == CL_SUCCESS, "error compiling programs:\n" + Tools::ErrorString( error_number));
    }

  } // namespace opencl
} // namespace bcl
