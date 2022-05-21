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
#include "model/bcl_model_approximator_support_vector_machine_multi_output.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_parameter_check_ranged.h"
#include "linal/bcl_linal_symmetric_eigensolver.h"
#include "model/bcl_model_kappa_nearest_neighbor.h"
#include "model/bcl_model_objective_function_constant.h"
#include "model/bcl_model_support_vector_machine_multi_output.h"
#include "util/bcl_util_stopwatch.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {
    const size_t ApproximatorSupportVectorMachineMultiOutput::s_Lower_Bound( 0); //!<  0 indicates lower boundary
    const size_t ApproximatorSupportVectorMachineMultiOutput::s_Upper_Bound( 1); //!<  1 indicates upper boundary
    const size_t ApproximatorSupportVectorMachineMultiOutput::s_Free( 2);        //!<  2 indicates no boundary
    const float ApproximatorSupportVectorMachineMultiOutput::m_EPS_A( 1e-12);

    // add the interface to the set of known implementations
    const util::SiPtr< const util::ObjectInterface> ApproximatorSupportVectorMachineMultiOutput::s_Instance
    (
      util::Enumerated< ApproximatorBase>::AddInstance( new ApproximatorSupportVectorMachineMultiOutput())
    );

    //! @brief default constructor
    ApproximatorSupportVectorMachineMultiOutput::ApproximatorSupportVectorMachineMultiOutput() :
      m_CostParameterC( 0.0),
      m_Status(),
      m_Alpha(),
      m_Gradient(),
      m_Bias(),
      m_Model( new SupportVectorMachineMultiOutput),
      m_ObjectiveFunction( GetDefaultObjectiveFunction()),
      m_NumberIterations( 0),
      m_NumberCurrentSupportVectors( 0),
      m_OptimizationGapThreshold( 0.1),
      m_OptimizationGap( 1.0),
      m_LastObjectiveFunctionValue( 0.0),
      m_Kappa( 1),
      m_CoordinateSystems(),
      m_CoordinateSystemsTransposed(),
      m_Epsilon( 0.001),
      m_SelfKernelValues()
    {
      ApproximatorBase::SetObjectiveFunction( GetDefaultObjectiveFunction());
    }

    //! @brief Iterate for Sequential Minimal Optimization Learning Algorithm
    //! @param COST_PARAMETER_C penalty parameter c for svm regression training
    //! @param MODEL initial support vector model
    //! @param TRAINING_DATA training data set of choice
    //! @param NUMBER_ITERATIONS
    ApproximatorSupportVectorMachineMultiOutput::ApproximatorSupportVectorMachineMultiOutput
    (
      const float COST_PARAMETER_C,
      const util::ShPtr< SupportVectorMachineMultiOutput> &MODEL,
      util::ShPtr< descriptor::Dataset> &TRAINING_DATA,
      const size_t NUMBER_ITERATIONS
    ) :
      m_CostParameterC( COST_PARAMETER_C),
      m_Status(),
      m_Alpha(),
      m_Gradient(),
      m_Bias(),
      m_Model( MODEL.IsDefined() ? MODEL : util::ShPtr< SupportVectorMachineMultiOutput>( new SupportVectorMachineMultiOutput())),
      m_ObjectiveFunction( GetDefaultObjectiveFunction()),
      m_NumberIterations( NUMBER_ITERATIONS),
      m_NumberCurrentSupportVectors( 0),
      m_OptimizationGapThreshold( 0.1),
      m_OptimizationGap( 1.0),
      m_LastObjectiveFunctionValue( 0.0),
      m_Kappa( 1),
      m_CoordinateSystems(),
      m_CoordinateSystemsTransposed(),
      m_Epsilon( 0.001),
      m_SelfKernelValues()
    {
      ApproximatorBase::SetObjectiveFunction( GetDefaultObjectiveFunction());

      // only used in conjunction with IteratateFromFile
      SetTrainingContinued( false);

      // set and rescale training data set
      SetTrainingData( TRAINING_DATA);
    }

    //! @brief copy constructor
    //! @return a new ApproximatorSupportVectorMachineMultiOutput copied from this instance
    ApproximatorSupportVectorMachineMultiOutput *ApproximatorSupportVectorMachineMultiOutput::Clone() const
    {
      return new ApproximatorSupportVectorMachineMultiOutput( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ApproximatorSupportVectorMachineMultiOutput::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &ApproximatorSupportVectorMachineMultiOutput::GetAlias() const
    {
      static const std::string s_Name( "MultiOutputSVM");
      return s_Name;
    }

    //! @brief set training data set for a specific iterate in approximater framework
    //! @param DATA training data set
    void ApproximatorSupportVectorMachineMultiOutput::SetTrainingData( util::ShPtr< descriptor::Dataset> &DATA)
    {
      // if the number of possible support vectors added in every iteration exceeds the number of training data
      // note: in every iteration can be 2 support vectors added
      if( 2 * m_NumberIterations > DATA->GetSize())
      {
        // set number of internal iterations to 1
        m_NumberIterations = 1;

        BCL_MessageDbg
        (
          "Reset number of internal iterations since one iterations adds "
          "more support vectors then training data available"
        );
      }

      m_TrainingData = DATA;
      DATA->GetFeatures().Rescale( SupportVectorMachineMultiOutput::s_DefaultInputRange, RescaleFeatureDataSet::e_AveStd);
      DATA->GetResults().Rescale( SupportVectorMachineMultiOutput::s_DefaultInputRange, RescaleFeatureDataSet::e_AveStd);

      // set rescale functions and kernel svm model used in iterative training process
      m_Model = util::ShPtr< SupportVectorMachineMultiOutput>
      (
        new SupportVectorMachineMultiOutput( m_Model->GetKernel(), GetRescaleFeatureDataSet(), GetRescaleResultDataSet())
      );

      // assemble SV Model
      m_Bias = linal::Vector< float>( m_TrainingData->GetResultSize(), float( 0.0));
      m_Model->SetBias( m_Bias);

      BCL_MessageStd( "before InitializeMemberVectorsForTraining");
      // if the current instance is NOT used through IterateInterfaceFromFile
      // then initialize the Interate properly
      if( !IsTrainingContinued())
      {
        BCL_MessageStd( "InitializeMemberVectorsForTraining");
        // initialize model and LLT coordinate systems
        InitializeMemberVectorsForTraining();
        ComputeCoordinateSystems();
      }
    }

    //! @brief construct a model from the current iterate
    //! @return shptr to the new model interface
    util::ShPtr< Interface> ApproximatorSupportVectorMachineMultiOutput::GetCurrentModel() const
    {
      return m_Model;
    }

    //! @brief returns the current approximation
    //! @return current argument result pair
    const util::ShPtr< storage::Pair< util::ShPtr< Interface>, float> >
      ApproximatorSupportVectorMachineMultiOutput::GetCurrentApproximation() const
    {
      util::ShPtr< storage::Pair< util::ShPtr< Interface>, float> > current_model;

      if
      (
        !m_ObjectiveFunction->GetImplementation().IsDefined()
        || m_ObjectiveFunction->GetImplementation().GetAlias() == "Constant"
      )
      {
        // create final pair with model and objective function evaluation
        current_model = util::ShPtr< storage::Pair< util::ShPtr< Interface>, float> >
        (
          new storage::Pair< util::ShPtr< Interface>, float>( m_Model.HardCopy(), m_OptimizationGap)
        );
      }
      else
      {
        // create final pair with model and objective function evaluation
        current_model = util::ShPtr< storage::Pair< util::ShPtr< Interface>, float> >
        (
          new storage::Pair< util::ShPtr< Interface>, float>( m_Model.HardCopy(), m_LastObjectiveFunctionValue)
        );
      }

      return current_model;
    }

    //! @brief conducts the next approximation step and stores the approximation
    void ApproximatorSupportVectorMachineMultiOutput::Next()
    {
      util::Stopwatch timer;

      if
      (
        m_OptimizationGap > m_OptimizationGapThreshold
        && m_NumberCurrentSupportVectors < float( 0.9) * m_TrainingData->GetSize()
      )
      {
        // iterate for a number of internal iterations
        for( size_t counter( 0); counter < m_NumberIterations; ++counter)
        {
           IterationStep();
        }

        // postprocess model and determine final support vectors
        m_OptimizationGap = FinalizeSupportVectorModel();

        m_NumberCurrentSupportVectors = m_Model->GetNumberSupportVectors();
        m_LastObjectiveFunctionValue = m_ObjectiveFunction->operator ()( util::ShPtr< Interface>( m_Model));
      }

      this->Track( GetCurrentApproximation());

      BCL_MessageStd
      (
        + " #SV: " + util::Format()( m_NumberCurrentSupportVectors)
        + " gap: " + util::Format()( m_OptimizationGap)
        + " threshold: " + util::Format()( m_OptimizationGapThreshold)
        + " time: " + util::Format()( timer.GetProcessDuration().GetTimeAsHourMinuteSecondMilliSeconds())
      );
    }

    //! @brief iterates one cycle and returns ShPtr to pair of resultant argument and corresponding score
    float ApproximatorSupportVectorMachineMultiOutput::IterationStep()
    {
      size_t first_vector_index( 0);
      size_t second_vector_index( 0);

      // determine index i and j for the two data vectors for which the quadratic problem has to be solved
      DetermineFeatureVectorCombination( first_vector_index, second_vector_index);

      // Solving Quadratic Problem sub problems for two given data vectors
      return SolveQuadraticProblemSubProblem( first_vector_index, second_vector_index, true);
    }

    //! @brief initialize a default objective function constant as default objective function
    util::ShPtr< ObjectiveFunctionWrapper> ApproximatorSupportVectorMachineMultiOutput::GetDefaultObjectiveFunction()
    {
      return util::ShPtr< ObjectiveFunctionWrapper>
      (
        new ObjectiveFunctionWrapper
        (
          util::Implementation< ObjectiveFunctionInterface>
          (
            new ObjectiveFunctionConstant
            (
              std::numeric_limits< float>::max(),
              opti::e_SmallerEqualIsBetter
            )
          )
        )
      );
    }

    //! @brief evaluates whether the approximation can continue
    //! @return true, if the approximation can continue - otherwise false
    bool ApproximatorSupportVectorMachineMultiOutput::CanContinue() const
    {
      return
        m_OptimizationGap > m_OptimizationGapThreshold
        && m_NumberCurrentSupportVectors < float( 0.9) * m_TrainingData->GetSize();
    }

    //! @brief finalize support vector model and determine final SVs, Alphas and Bias
    float ApproximatorSupportVectorMachineMultiOutput::FinalizeSupportVectorModel()
    {
      // counter of support vectors
      size_t beta_count( 0);

      // counter for iteration
      size_t current_vector_index( 0);

      //calculate bias
      const linal::Vector< float> bias( CalculateBias());
      m_Bias = bias;

      // determine optimization gap = sqrt( sum( (error - bias)^2) / ( number_outputs * number_features))
      // this is similar to standard derivation and measures the deviation of predictions and outputs
      // sum of error vectors squares
      float error_sum( 0.0);

      // go over all error/gradient vectors and square them using dot product
      for
      (
          storage::Vector< linal::Vector< float> >::const_iterator itr_error( m_Gradient.Begin()),
          itr_error_end( m_Gradient.End());
          itr_error != itr_error_end;
          ++itr_error
      )
      {
        const linal::Vector< float> reduced_gradient( *itr_error - bias);
        error_sum += reduced_gradient * reduced_gradient;
      }

      // to get a standard-derivation-like value divide by the number of outputs and take the sqrt
      // it is independent of the number of outputs
      const float optimization_gap
      (
        math::Sqrt( error_sum / float( m_TrainingData->GetSize() * m_TrainingData->GetResultSize()))
      );

      // count alpha vector elements that hit C boundary
      // there should be a reasonable number of bound alphas if there are too much or too less adjust m_CostParameterC
      size_t bound_alphas( 0);
      for( size_t index_alpha( 0); index_alpha < m_TrainingData->GetFeatureSize(); ++index_alpha)
      {
        for( size_t alpha_element( 0); alpha_element < m_TrainingData->GetResultSize(); ++alpha_element)
        {
          if( !IsFree( index_alpha, alpha_element))
          {
            ++bound_alphas;
          }
        }
      }

      BCL_Message( util::Message::e_Standard, "Bound Alphas: " + util::Format()( bound_alphas));

      // final vector with all betas of support vectors. beta = transposed coordinate system * alpha
      storage::Vector< linal::Vector< float> > betas;

      // list of indizes of the support vectors
      storage::List< size_t> sv_indices;

      // assembly of the final beta vector of the lagrange multipliers
      storage::Vector< linal::Matrix< float> >::const_iterator iter_coord_sys( m_CoordinateSystems.Begin());

      // go over all alpha vectors
      for
      (
        storage::Vector< linal::Vector< float> >::const_iterator iter_alpha( m_Alpha.Begin()),
        iter_alpha_end( m_Alpha.End());
        iter_alpha != iter_alpha_end;
        ++iter_alpha, ++iter_coord_sys
      )
      {
        // iterate over all alpha elements and set them to zero if their value is smaller than the threshold m_Epsilon
        linal::Vector< float> current_alpha( *iter_alpha);
        for
        (
            linal::Vector< float>::iterator iter_alpha_element( current_alpha.Begin()),
            iter_alpha_element_end( current_alpha.End());
            iter_alpha_element != iter_alpha_element_end;
            ++iter_alpha_element
        )
        {
          // check whether alpha element is in epsilon tube
          if( *iter_alpha_element < m_Epsilon && *iter_alpha_element > -m_Epsilon)
          {
            *iter_alpha_element = 0;
          }
        }
        // compute the beta vector
        const linal::Vector< float> current_beta( *iter_coord_sys * ( current_alpha));

        // make it support vector if at least one element is nonzero
        bool is_support_vector( false);

        // iterate through all betas
        for
        (
            linal::Vector< float>::const_iterator iter_beta( current_beta.Begin()),
            iter_beta_end( current_beta.End());
            iter_beta != iter_beta_end;
            ++iter_beta
        )
        {
          // if beta is not zero then correspondent vector is a support vector
          if( *iter_beta > m_EPS_A || *iter_beta < -m_EPS_A)
          {
            is_support_vector = true;
          }
        }
        // store beta of support vector in variable betas
        if( is_support_vector)
        {
          betas.PushBack( current_beta);
          sv_indices.PushBack( current_vector_index);
          ++beta_count;
        }

        // increase iteration counter
        ++current_vector_index;
      }

      // matrix containing support vectors
      linal::Matrix< float> sv_matrix
      (
        sv_indices.GetSize(),                             // rows
        m_TrainingData->GetFeatureSize(),                 // cols
        float( 0)                                         // default value
      );

      // reset counter
      current_vector_index = 0;

      // fill support vector matrix with vectors by support vector index
      for
      (
        storage::List< size_t>::const_iterator itr_sv( sv_indices.Begin()), itr_sv_end( sv_indices.End());
        itr_sv != itr_sv_end;
        ++itr_sv, ++current_vector_index
      )
      {
        // get current support vector
        const FeatureReference< float> &support_vector( m_TrainingData->GetFeatures()( *itr_sv));
        // copy support vector to support vector matrix
        std::copy( support_vector.Begin(), support_vector.End(), sv_matrix[ current_vector_index]);
      }

      // assemble SV Model
      m_Model->SetBeta( betas);
      m_Model->SetBias( bias);
      m_Model->SetSupportVectors( FeatureDataSet< float>( sv_matrix));
      m_Model->SetNumberSupportVectors( m_Model->GetSupportVectors().GetNumberFeatures());

      return optimization_gap;
    }

    //! @brief update status of one complete alpha vector according to its alpha vector
    //! @param FEATURE_VECTOR_I a position of a feature vector in m_Status
    void ApproximatorSupportVectorMachineMultiOutput::UpdateAlphaStatus( const size_t &FEATURE_VECTOR_I)
    {
      const linal::Vector< float> &alpha_i( m_Alpha( FEATURE_VECTOR_I));
      storage::Vector< size_t> &status_i( m_Status( FEATURE_VECTOR_I));
      linal::Vector< float>::const_iterator itr_alpha( alpha_i.Begin());
      storage::Vector< size_t>::iterator itr_stat( status_i.Begin()), itr_stat_end( status_i.End());

      // update boundary status for every element of the alpha vector
      for( ; itr_stat != itr_stat_end; ++itr_stat, ++itr_alpha)
      {
        if( *itr_alpha >= m_CostParameterC - m_EPS_A)
        {
          *itr_stat = s_Upper_Bound;
        }
        else if( *itr_alpha <= m_EPS_A - m_CostParameterC)
        {
          *itr_stat = s_Lower_Bound;
        }
        else
        {
          *itr_stat = s_Free;
        }
      }
    }

    //! read NeuralNetwork from std::istream
    std::istream &ApproximatorSupportVectorMachineMultiOutput::Read( std::istream &ISTREAM)
    {
      storage::Vector< float>::s_Instance.IsDefined();

      // read members
      io::Serialize::Read( m_CostParameterC, ISTREAM);
      io::Serialize::Read( m_Status, ISTREAM);
      io::Serialize::Read( m_Alpha, ISTREAM);
      io::Serialize::Read( m_Gradient, ISTREAM);
      io::Serialize::Read( m_Bias, ISTREAM);
      io::Serialize::Read( m_Model, ISTREAM);
      io::Serialize::Read( m_NumberIterations, ISTREAM);
      io::Serialize::Read( m_OptimizationGapThreshold, ISTREAM);
      io::Serialize::Read( m_OptimizationGap, ISTREAM);
      io::Serialize::Read( m_LastObjectiveFunctionValue, ISTREAM);
      io::Serialize::Read( m_Kappa, ISTREAM);
      io::Serialize::Read( m_CoordinateSystems, ISTREAM);
      io::Serialize::Read( m_CoordinateSystemsTransposed, ISTREAM);
      io::Serialize::Read( m_Epsilon, ISTREAM);
      io::Serialize::Read( m_SelfKernelValues, ISTREAM);

      // return
      return ISTREAM;
    }

    //! write NeuralNetwork into std::ostream
    std::ostream &ApproximatorSupportVectorMachineMultiOutput::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_CostParameterC, OSTREAM, INDENT) << std::endl;
      io::Serialize::Write( m_Status, OSTREAM, INDENT) << std::endl;
      io::Serialize::Write( m_Alpha, OSTREAM, INDENT) << std::endl;
      io::Serialize::Write( m_Gradient, OSTREAM, INDENT) << std::endl;
      io::Serialize::Write( m_Bias, OSTREAM, INDENT) << std::endl;
      io::Serialize::Write( m_Model, OSTREAM, INDENT) << std::endl;
      io::Serialize::Write( m_NumberIterations, OSTREAM, INDENT) << std::endl;
      io::Serialize::Write( m_OptimizationGapThreshold, OSTREAM, INDENT) << std::endl;
      io::Serialize::Write( m_OptimizationGap, OSTREAM, INDENT) << std::endl;
      io::Serialize::Write( m_LastObjectiveFunctionValue, OSTREAM, INDENT) << std::endl;
      io::Serialize::Write( m_Kappa, OSTREAM, INDENT) << std::endl;
      io::Serialize::Write( m_CoordinateSystems, OSTREAM, INDENT) << std::endl;
      io::Serialize::Write( m_CoordinateSystemsTransposed, OSTREAM, INDENT) << std::endl;
      io::Serialize::Write( m_Epsilon, OSTREAM, INDENT) << std::endl;
      io::Serialize::Write( m_SelfKernelValues, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

    //! @brief method checks whether a feature_vector i reached the upper boundary
    //! @param FEATURE_VECTOR_I a position of a feature vector in m_Status
    //! @param ELEMENT position of the element - corresponds to the output column
    //! @return const bool that indicates whether FEATURE_VECTOR_I reached upper boundary
    inline bool ApproximatorSupportVectorMachineMultiOutput::IsUpperBound( const size_t &FEATURE_VECTOR_I, const size_t &ELEMENT) const
    {
      return m_Status( FEATURE_VECTOR_I)( ELEMENT) == s_Upper_Bound;
    }
    //! @brief checks whether a feature_vector i reached the lower boundary
    //! @param FEATURE_VECTOR_I a position of a feature vector in m_Status
    //! @param ELEMENT position of the element - corresponds to the output column
    //! @return const bool that indicates whether FEATURE_VECTOR_I reached lower boundary
    inline bool ApproximatorSupportVectorMachineMultiOutput::IsLowerBound( const size_t &FEATURE_VECTOR_I, const size_t &ELEMENT) const
    {
      return m_Status( FEATURE_VECTOR_I)( ELEMENT) == s_Lower_Bound;
    }

    //! @brief checks whether a feature_vector i is not bound and between the boundaries
    //! @param FEATURE_VECTOR_I a position of a feature vector in m_Status
    //! @param ELEMENT position of the element - corresponds to the output column
    //! @return const bool that indicates whether FEATURE_VECTOR_I reached no boundary
    inline bool ApproximatorSupportVectorMachineMultiOutput::IsFree( const size_t &FEATURE_VECTOR_I, const size_t &ELEMENT) const
    {
      return m_Status( FEATURE_VECTOR_I)( ELEMENT) == s_Free;
    }

    //! @brief Initialize member vectors and variables to prepare SVR training
    void ApproximatorSupportVectorMachineMultiOutput::InitializeMemberVectorsForTraining()
    {
      // Initialization of variables and vectors
      const size_t number_feature_vectors( m_TrainingData->GetSize());
      const size_t result_size( m_TrainingData->GetResultSize());

      // Initialize lagrange multiplier vectors
      m_Alpha = storage::Vector< linal::Vector< float> >
      (
        number_feature_vectors,
        linal::Vector< float>( result_size, float( 0.0))
      );

      //initialize alpha status
      m_Status = storage::Vector< storage::Vector< size_t> >
      (
        number_feature_vectors, storage::Vector< size_t>( result_size, size_t( 0))
      );

      //initialize alpha status
      for( size_t alpha_vector_index( 0); alpha_vector_index < number_feature_vectors; ++alpha_vector_index)
      {
        UpdateAlphaStatus( alpha_vector_index);
      }

      //compute the kernel value of each input vector with itself and store them in m_SelfKernelValues
      float kernel_value( 0.0);
      size_t progress( 0);
      m_SelfKernelValues = linal::Vector< float>( number_feature_vectors, float( 0.0));
      for
      (
        linal::Vector< float>::iterator iter_begin( m_SelfKernelValues.Begin()),
        iter_end( m_SelfKernelValues.End());
        iter_begin != iter_end;
        ++iter_begin, ++progress
      )
      {
        // get current feature
        const FeatureReference< float> &item_ref_at_position( m_TrainingData->GetFeatures()( progress));

        // compute self kernel value
        kernel_value = m_Model->GetKernel()->operator()( item_ref_at_position, item_ref_at_position);

        *iter_begin = kernel_value;
      }

      //initialize Error Vectors
      //Because all alphas are 0 in the beginning, m_Gradient = result
      //bias doesn't matter during training and gradients must not contain bias for bias calculation
      m_Gradient.AllocateMemory( number_feature_vectors);

      for( size_t i( 0); i < number_feature_vectors; ++i)
      {
        m_Gradient.PushBack( m_TrainingData->GetResults()( i));
      }

      BCL_Message( util::Message::e_Debug, " Initialization.. complete");
    }

    //! @brief determine the bias value for Support Vector Regression Model
    //! @return calculated bias value
    linal::Vector< float> ApproximatorSupportVectorMachineMultiOutput::CalculateBias() const
    {
      // initialization
      const size_t &result_size( m_TrainingData->GetResultSize());
      linal::Vector< float> bias( result_size, float( 0.0));
      int nr_free_gradients( 0);
      linal::Vector< float> sum_free_gradients( result_size, float( 0.0));

      // initialization of progress counter
      size_t progress( 0);

      // iterate over all gradients and sum gradients of non bound alphas
      for
      (
        storage::Vector< linal::Vector< float> >::const_iterator iter_gradient( m_Gradient.Begin()),
        iter_end_gradient( m_Gradient.End());
        iter_gradient != iter_end_gradient;
        ++iter_gradient, ++progress
      )
      {
        bool is_free( true);

        for( size_t element_id( 0); element_id < result_size; ++element_id)
        {
          // check for boundary conditions and update upper or lower boundaries accordingly
          if( !IsFree( progress, element_id))
          {
            is_free = false;
          }
        }
        if( is_free)
        {
          ++nr_free_gradients;
          sum_free_gradients += *iter_gradient;
        }
      }

      BCL_Assert( nr_free_gradients > 0, "All alpha vectors hit bounds. Adjust Parameter C");

      // calculate bias
      // iterate over all outputs
      for
      (
          linal::Vector< float>::iterator iter_bias( bias.Begin()),
          iter_sum_free( sum_free_gradients.Begin()),
          iter_end_bias( bias.End());
          iter_bias != iter_end_bias;
          ++iter_bias, ++iter_sum_free
      )
      {
        //if there are free elements take the average bias
        *iter_bias = *iter_sum_free / float( nr_free_gradients);
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
    void ApproximatorSupportVectorMachineMultiOutput::DetermineFeatureVectorCombination
    (
      size_t &FIRST_VECTOR_INDEX,
      size_t &SECOND_VECTOR_INDEX
    )
    {
      // initialize stopwatch
      util::Stopwatch timer;
      timer.Reset();
      timer.Start();

      // get current feature
      const float &number_features( m_TrainingData->GetFeaturesPtr()->GetNumberFeatures());

      // initialize
      float gradient_vector_i_sqr( -1);
      size_t feature_index_i( 0);
      bool found_feature_i( false);
      size_t process( 0);

      // determine i such that:
      // error vector transformed in local coordinate system has a maximum value,
      // bound hitting elements are ignored
      for
      (
        storage::Vector< linal::Vector< float> >::const_iterator
          iter_gradient( m_Gradient.Begin()),
          iter_end_gradient( m_Gradient.End());
        iter_gradient != iter_end_gradient;
        ++iter_gradient, ++process
      )
      {
        // initialize
        size_t element( 0);
        float sqr_sum( 0.0);

        // transform into local coordinate system
        const linal::Vector< float> trans_grad( m_CoordinateSystemsTransposed( process) * ( *iter_gradient - m_Bias));

        //iterate over all elements of transformed gradient
        for
        (
            linal::Vector< float>::const_iterator
              iter_trans_grad_element( trans_grad.Begin()),
              iter_end_trans_grad_element( trans_grad.End());
            iter_trans_grad_element != iter_end_trans_grad_element;
            ++iter_trans_grad_element, ++element
        )
        {
          // if the alpha would be adjusted in a direction where the boundary is hit, this element doesn't count
          if( *iter_trans_grad_element > 0)
          {
            if( !IsUpperBound( process, element))
            {
              sqr_sum += math::Sqr( *iter_trans_grad_element);
            }
          }
          else
          {
            if( !IsLowerBound( process, element))
            {
              sqr_sum += math::Sqr( *iter_trans_grad_element);
            }
          }
        }

        // save the maximum transformed error vector
        if( sqr_sum > gradient_vector_i_sqr)
        {
          gradient_vector_i_sqr = sqr_sum;
          feature_index_i = process;
          found_feature_i = true;
        }
      }

      // this happens if all gradients are zero or all alphas hit bounds
      BCL_Assert( found_feature_i, "no Vector i found");

      //for choosing j we solve the problem for every vector and take the one with the largest change in
      //beta = trans_coord_system * alpha
      float best_delta_bata_value( 0);
      size_t best_j( 0);
      for( size_t j( 0); j < number_features; ++j)
      {
        //finds the solution but doesn't apply it
        float delta_beta( SolveQuadraticProblemSubProblem( feature_index_i, j, false));

        if( best_delta_bata_value < delta_beta)
        {
          best_delta_bata_value = delta_beta;
          best_j = j;
        }
      }

      const size_t index_j( best_j);

      FIRST_VECTOR_INDEX  = feature_index_i;
      SECOND_VECTOR_INDEX = index_j;

      BCL_Message
      (
        util::Message::e_Debug, "selected i:" + util::Format()( feature_index_i)
        + " selected j:" + util::Format()( index_j)
      );

      timer.Stop();
      const util::Time total_duration_so_far( timer.GetTotalTime());

      BCL_Message( util::Message::e_Debug, "Select Duration: " + util::Format()( total_duration_so_far));
    }

    //! @brief solve sub problem of quadratic problem according two feature vectors
    //! @param SVR_MODEL support vector machine model of interest
    //! @param FEATURE_VECTOR_I index of feature vector in MATRIX
    //! @param FEATURE_VECTOR_J index of feature vector in MATRIX
    //! @return a pair of vector of float with computed values of kernel function
    float ApproximatorSupportVectorMachineMultiOutput::SolveQuadraticProblemSubProblem
    (
      const size_t &FEATURE_VECTOR_I,
      const size_t &FEATURE_VECTOR_J,
      const bool &APPLY_SOLUTION
    )
    {
      //get result size
      const size_t result_size( m_TrainingData->GetResultSize());

      // compute kernel value of chosen feature vectors
      const float k_ij
      (
        m_Model->GetKernel()->operator()
        (
          m_TrainingData->GetFeaturesPtr()->operator ()( FEATURE_VECTOR_I),
          m_TrainingData->GetFeaturesPtr()->operator ()( FEATURE_VECTOR_J)
        )
      );

      // compute quad_coef eta
      float quad_coef( m_SelfKernelValues( FEATURE_VECTOR_I) - 2 * k_ij + m_SelfKernelValues( FEATURE_VECTOR_J));
      if( quad_coef <= m_EPS_A) //TDO use max()
      {
        quad_coef = m_EPS_A;
      }

      //get local coordinate_systems and the transposed matrices
      const linal::Matrix< float> &coord_sys_i( m_CoordinateSystems( FEATURE_VECTOR_I));
      const linal::Matrix< float> &coord_sys_j( m_CoordinateSystems( FEATURE_VECTOR_J));
      const linal::Matrix< float> &coord_sys_i_trans( m_CoordinateSystemsTransposed( FEATURE_VECTOR_I));
      const linal::Matrix< float> &coord_sys_j_trans( m_CoordinateSystemsTransposed( FEATURE_VECTOR_J));

      //get error vectors at i and j
      linal::Vector< float> error_i( m_Gradient( FEATURE_VECTOR_I));
      linal::Vector< float> error_j( m_Gradient( FEATURE_VECTOR_J));

      //get Alphas
      const linal::Vector< float> &alpha_i_old( m_Alpha( FEATURE_VECTOR_I));
      const linal::Vector< float> &alpha_j_old( m_Alpha( FEATURE_VECTOR_J));

      //get gamma used in the summation constraint
      const linal::Vector< float> gamma( coord_sys_i * alpha_i_old + coord_sys_j * alpha_j_old);

      //compute new unclipped alphas
      linal::Vector< float> alpha_i_new( result_size, float( 0.0));
      linal::Vector< float> alpha_j_new( result_size, float( 0.0));

      // compute the unclipped alpha_i
      alpha_i_new = alpha_i_old + coord_sys_i_trans * ( error_i - error_j) / quad_coef;

      // fit alpha_i into boundary constraints
      for
      (
        linal::Vector< float>::iterator
        itr_alpha = alpha_i_new.Begin(),
        itr_alpha_end = alpha_i_new.End();
        itr_alpha != itr_alpha_end;
        ++itr_alpha
      )
      {
        if( *itr_alpha <= -m_CostParameterC - m_EPS_A)
        {
          *itr_alpha = -m_CostParameterC;
        }
        if( *itr_alpha >= m_CostParameterC + m_EPS_A)
        {
          *itr_alpha = m_CostParameterC;
        }
      }

      // check for boundary violation and loop-ajust the alphas
      bool bounds_correct( false);

      // set a constant value for max iterations to fit alphas into bounds
      const size_t max_bound_iterations( 5);

      for( size_t bound_iterations( 0); bound_iterations < max_bound_iterations && !bounds_correct; ++bound_iterations)
      {
        // calculate alpha_j using the summation constraint
        alpha_j_new = coord_sys_j_trans * ( gamma - coord_sys_i * alpha_i_new);

        // check bounds of alpha_j
        bounds_correct = true;
        for
        (
            linal::Vector< float>::iterator
            itr = alpha_j_new.Begin(),
            itr_end = alpha_j_new.End();
            itr != itr_end;
            ++itr
        )
        {
          // alpha can be numerical unstable, that's why violations smaller than m_EPS_A are not taken care of
          if( *itr > m_CostParameterC + m_EPS_A)
          {
            *itr = m_CostParameterC;
            bounds_correct = false;
          }
          if( *itr < -m_CostParameterC - m_EPS_A)
          {
            *itr = -m_CostParameterC;
            bounds_correct = false;
          }
        }
        if( bounds_correct)
        {
          break;
        }
        // compute cliped alpha_i if bounds were violated
        alpha_i_new = coord_sys_i_trans * ( gamma - coord_sys_j * alpha_j_new);

        // check bounds of alpha_i
        bounds_correct = true;
        for
        (
            linal::Vector< float>::iterator
            itr = alpha_i_new.Begin(),
            itr_end = alpha_i_new.End();
            itr != itr_end;
            ++itr
        )
        {
          if( *itr > m_CostParameterC + m_EPS_A)
          {
            *itr = m_CostParameterC;
            bounds_correct = false;
          }
          if( *itr < -m_CostParameterC - m_EPS_A)
          {
            *itr = -m_CostParameterC;
            bounds_correct = false;
          }
        }
      }

      if( !bounds_correct)
      {
        BCL_Message( util::Message::e_Debug, "No correct bounds during BoundInterations");
        //fix bounds ignoring summation constraint
        for
        (
            linal::Vector< float>::iterator
            itr_i = alpha_i_new.Begin(),
            itr_j = alpha_j_new.Begin(),
            itr_i_end = alpha_i_new.End();
            itr_i != itr_i_end;
            ++itr_i, ++itr_j
        )
        {
          if( *itr_i > m_CostParameterC + m_EPS_A)
          {
            *itr_i = m_CostParameterC;
          }
          if( *itr_i < -m_CostParameterC - m_EPS_A)
          {
            *itr_i = -m_CostParameterC;
          }
          if( *itr_j > m_CostParameterC + m_EPS_A)
          {
            *itr_j = m_CostParameterC;
          }
          if( *itr_j < -m_CostParameterC - m_EPS_A)
          {
            *itr_j = -m_CostParameterC;
          }
        }
      }

      //calculate delta beta
      const linal::Vector< float> delta_beta( coord_sys_i * ( alpha_i_new - alpha_i_old));
      const float delta_beta_value( delta_beta * delta_beta);

      if( APPLY_SOLUTION)
      {
        // using Kernel Q_i(x) = K(bound_iterations,x)
        const linal::Vector< float> Q_i
        (
          m_Model->GetKernel()->GetInputVectorIKernelMatrixMultiOutput( *( m_TrainingData->GetFeaturesPtr()), FEATURE_VECTOR_I)
        );
        // using Kernel Q_j(x) = K(j,x)
        const linal::Vector< float> Q_j
        (
          m_Model->GetKernel()->GetInputVectorIKernelMatrixMultiOutput( *( m_TrainingData->GetFeaturesPtr()), FEATURE_VECTOR_J)
        );

        //update Error vectors
        //transform the delta_alpha into the local coordinate system
        const linal::Vector< float> delta_alpha_i( coord_sys_i * ( alpha_i_new - alpha_i_old));
        const linal::Vector< float> delta_alpha_j( coord_sys_j * ( alpha_j_new - alpha_j_old));
        linal::Vector< float>::const_iterator iter_begin_Qi( Q_i.Begin()), iter_begin_Qj( Q_j.Begin());

        for
        (
          storage::Vector< linal::Vector< float> >::iterator
            iter_begin_gradient( m_Gradient.Begin()),
            iter_end_gradient( m_Gradient.End());
          iter_begin_gradient != iter_end_gradient;
          ++iter_begin_gradient, ++iter_begin_Qi, ++iter_begin_Qj
        )
        {
          *iter_begin_gradient -= *iter_begin_Qi * delta_alpha_i + *iter_begin_Qj * delta_alpha_j;
        }
        //submit new Alphas to m_Alpha
        m_Alpha( FEATURE_VECTOR_I) = alpha_i_new;
        m_Alpha( FEATURE_VECTOR_J) = alpha_j_new;

        //update Alpha status
        UpdateAlphaStatus( FEATURE_VECTOR_I);
        UpdateAlphaStatus( FEATURE_VECTOR_J);
      }

      return delta_beta_value;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer ApproximatorSupportVectorMachineMultiOutput::GetSerializer() const
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
        GetDefaultObjectiveFunction()->GetImplementation().GetString()
      );
      parameters.AddInitializer
      (
        "kernel",
        "kernel used to map pairs of features onto a hyperplane",
        io::Serialization::GetAgent( &m_Model->GetKernel()),
        "RBF(gamma=0.5)"
      );

      parameters.AddInitializer
      (
        "iterations",
        "# of iterations used internally to improve the optimization gap",
        io::Serialization::GetAgent( &m_NumberIterations),
        "1"
      );
      parameters.AddInitializer
      (
        "cost",
        "controls the trade off between allowing training errors and forcing rigid margins; high values may lead to better training values, but risks overtraining",
        io::Serialization::GetAgent( &m_CostParameterC),
        "0.0"
      );
      parameters.AddInitializer
      (
        "gap_threshold",
        "optimization gap threshold - the training process will stop when threshold is reached, higher values lead to more training cycles ",
        io::Serialization::GetAgent( &m_OptimizationGapThreshold),
        "0.1"
      );

      parameters.AddInitializer
      (
        "kappa",
        "Number of nearest neighbors used for computation of locally linear transformation coordinate system",
        io::Serialization::GetAgent( &m_Kappa),
        "1"
      );

      parameters.AddInitializer
      (
        "epsilon",
        "Range within errors are ignored. The smaller the more accurate is the prediction and the more support vectors are needed ",
        io::Serialization::GetAgent( &m_Epsilon),
        "0.0001"
      );

      return parameters;
    }

    //! @brief computes for each output vector a local coordinate system used in the multi output method
    void ApproximatorSupportVectorMachineMultiOutput::ComputeCoordinateSystems()
    {
      const size_t result_size( m_TrainingData->GetResultSize());
      const size_t train_size( m_TrainingData->GetSize());

      BCL_Assert( m_Kappa < train_size, "Kappa is not smaller than Training set");
      // If kappa is one or less the resulting coordianate systems are set to the identity matrix
      if( m_Kappa < 2)
      {
        linal::Matrix< float> identity( result_size, result_size, float( 0));
        for( size_t i( 0); i < result_size; ++i)
        {
          identity( i, i) = 1;
        }
        for( size_t i( 0); i < train_size; ++i)
        {
          m_CoordinateSystems.PushBack( identity);
          m_CoordinateSystemsTransposed.PushBack( identity);
        }
      }
      //else compute them from kappa neighbors
      else
      {
        linal::Matrix< float> neighbor_distance( result_size, m_Kappa, float( 0.0));

        // for each result vector
        for( size_t current_res_vector = 0; current_res_vector < train_size; ++current_res_vector)
        {
          //find the number of equal vectors
          size_t number_equals( 0);
          for( size_t i( 0); i < train_size; ++i)
          {
            if
            (
              math::EqualWithinTolerance
              (
                m_TrainingData->GetResults()( i),
                m_TrainingData->GetResults()( current_res_vector)
              )
            )
            {
              ++number_equals;
            }
          }

          // if all outputs are equal
          BCL_Assert( number_equals < train_size, "No different Outputs");
          size_t clipped_kappa( m_Kappa);
          if( m_Kappa + number_equals > train_size)
          {
            BCL_Message
            (
              util::Message::e_Standard,
              "On vector " + util::Format()( current_res_vector) +
              " are less than kappa different vectors"
            );
            // clipp Kappa to half the different vectors, at least one
            clipped_kappa = ( train_size - number_equals) / 2;
          }

          // get the Kappa nearest neighbors
          storage::Vector< storage::Pair< float, size_t> > neighbors
          (
            KappaNearestNeighbor::FindWithoutRescaling
            (
              m_TrainingData->GetResults()( current_res_vector),
              clipped_kappa + number_equals, m_TrainingData->GetResults()
            )
          );

          //get average
          linal::Vector< float> average( result_size, float( 0.0));

          for( size_t j( 0); j < clipped_kappa; ++j)
          {
            average += linal::Vector< float>
            (
              result_size,
              m_TrainingData->GetResults().GetMatrix()[ neighbors( j + number_equals).Second()]
            ) / float( clipped_kappa);
          }

          // fill matrix with difference to the average
          for( size_t i( 0); i < result_size; ++i)
          {
            for( size_t j( 0); j < clipped_kappa; ++j)
            {
              neighbor_distance( i, j) = m_TrainingData->GetResults()( neighbors( j + number_equals).Second())( i) - average( i);
            }
          }
          //do singular value decomposition
          //this equals getting the eigenvectors of the covariance matrix
          // We need the left-singular vectors as a square matrix. That's equal to compute the right-singular vectors of the transposed matrix
          linal::Matrix< float> u( 1, 1), v, vt( 1, 1);
          linal::SingularValueDecomposition( neighbor_distance.Transposed(), vt, u);
          v = vt.Transposed();
          //copy V into m_CoordinateSystems
          m_CoordinateSystems.PushBack( v);

          //copy transposed matrix into m_CoordinateSystemsTransposed
          m_CoordinateSystemsTransposed.PushBack( vt);
        }
      }
      BCL_Message( util::Message::e_Debug, "CoordSys:" + util::Format()( m_CoordinateSystems));
    }

  } // namespace model
} // namespace bcl
