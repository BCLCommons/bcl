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
#include "model/bcl_model_approximator_support_vector_machine.h"

// includes from bcl - sorted alphabetically
#include "model/bcl_model_objective_function_constant.h"
#include "model/bcl_model_support_vector_machine.h"
#include "util/bcl_util_stopwatch.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {
    const size_t ApproximatorSupportVectorMachine::s_Lower_Bound( 0); //!<  0 indicates lower boundary
    const size_t ApproximatorSupportVectorMachine::s_Upper_Bound( 1); //!<  1 indicates upper boundary
    const size_t ApproximatorSupportVectorMachine::s_Free( 2);        //!<  2 indicates no boundary
    const float ApproximatorSupportVectorMachine::m_EPS_A( 1e-12);
    const float ApproximatorSupportVectorMachine::m_P( 0.1);

    // add the interface to the set of known implementations
    const util::SiPtr< const util::ObjectInterface> ApproximatorSupportVectorMachine::s_Instance
    (
      util::Enumerated< ApproximatorBase>::AddInstance( new ApproximatorSupportVectorMachine())
    );

    //! @brief default constructor
    ApproximatorSupportVectorMachine::ApproximatorSupportVectorMachine() :
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
      m_Model( new SupportVectorMachine),
      m_NumberIterations( 0),
      m_NumberCurrentSupportVectors( 0),
      m_OptimizationGapThreshold( 0.1),
      m_OptimizationGap( 1.0),
      m_LastObjectiveFunctionValue( 0.0)
    {
      ApproximatorBase::SetObjectiveFunction( GetDefaultObjectiveFunction());
    }

    //! @brief Iterate for Sequential Minimal Optimization Learning Algorithm
    //! @param COST_PARAMETER_C penalty parameter c for svm regression training
    //! @param MODEL initial support vector model
    //! @param TRAINING_DATA training data set of choice
    //! @param NUMBER_ITERATIONS
    ApproximatorSupportVectorMachine::ApproximatorSupportVectorMachine
    (
      const float COST_PARAMETER_C,
      const util::ShPtr< SupportVectorMachine> &MODEL,
      util::ShPtr< descriptor::Dataset> &TRAINING_DATA,
      const size_t NUMBER_ITERATIONS
    ) :
      m_CostParameterC( COST_PARAMETER_C),
      m_Status( storage::Vector< size_t>()),
      m_Alpha( storage::Vector< float>()),
      m_Gradient( storage::Vector< float>()),
      m_GradientBar( storage::Vector< float>()),
      m_Bias( storage::Vector< float>()),
      m_Signs( storage::Vector< int>()),
      m_Labels( storage::Vector< float>()),
      m_ActiveSize( 0),
      m_ProbLength( 2 * TRAINING_DATA->GetSize()),
      m_Model( MODEL.IsDefined() ? MODEL : util::ShPtr< SupportVectorMachine>( new SupportVectorMachine())),
      m_NumberIterations( NUMBER_ITERATIONS),
      m_NumberCurrentSupportVectors( 0),
      m_OptimizationGapThreshold( 0.1),
      m_OptimizationGap( 1.0),
      m_LastObjectiveFunctionValue( 0.0)
    {
      ApproximatorBase::SetObjectiveFunction( GetDefaultObjectiveFunction());

      // only used in conjunction with IteratateFromFile
      SetTrainingContinued( false);

      // set and rescale training data set
      SetTrainingData( TRAINING_DATA);
    }

    //! @brief copy constructor
    //! @return a new ApproximatorSupportVectorMachine copied from this instance
    ApproximatorSupportVectorMachine *ApproximatorSupportVectorMachine::Clone() const
    {
      return new ApproximatorSupportVectorMachine( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ApproximatorSupportVectorMachine::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &ApproximatorSupportVectorMachine::GetAlias() const
    {
      static const std::string s_Name( "SupportVectorMachine");
      return s_Name;
    }

    //! @brief set training data set for a specific iterate in approximater framework
    //! @param DATA training data set
    void ApproximatorSupportVectorMachine::SetTrainingData( util::ShPtr< descriptor::Dataset> &DATA)
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
      DATA->GetFeatures().Rescale( SupportVectorMachine::s_DefaultInputRange, RescaleFeatureDataSet::e_AveStd);
      DATA->GetResults().Rescale( SupportVectorMachine::s_DefaultInputRange, RescaleFeatureDataSet::e_AveStd);

      // set rescale functions and kernel svm model used in iterative training process
      m_Model = util::ShPtr< SupportVectorMachine>
      (
        new SupportVectorMachine( m_Model->GetKernel(), GetRescaleFeatureDataSet(), GetRescaleResultDataSet())
      );

      BCL_MessageStd( "before InitializeMemberVectorsForTraining");
      // if the current instance is NOT used through IterateInterfaceFromFile
      // then initialize the Interate properly
      if( !IsTrainingContinued())
      {
        BCL_MessageStd( "InitializeMemberVectorsForTraining");
        // initialize model
        InitializeMemberVectorsForTraining();
      }
    }

    //! @brief construct a model from the current iterate
    //! @return shptr to the new model interface
    util::ShPtr< Interface> ApproximatorSupportVectorMachine::GetCurrentModel() const
    {
      return m_Model;
    }

    //! @brief returns the current approximation
    //! @return current argument result pair
    const util::ShPtr< storage::Pair< util::ShPtr< Interface>, float> >
      ApproximatorSupportVectorMachine::GetCurrentApproximation() const
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
    void ApproximatorSupportVectorMachine::Next()
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
          m_OptimizationGap = IterationStep();
        }

        // postprocess model and determine final support vectors
        FinalizeSupportVectorModel();

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
    float ApproximatorSupportVectorMachine::IterationStep()
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

    //! @brief initialize a default objective function constant as default objective function
    util::ShPtr< ObjectiveFunctionWrapper> ApproximatorSupportVectorMachine::GetDefaultObjectiveFunction()
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
    bool ApproximatorSupportVectorMachine::CanContinue() const
    {
      return
        m_OptimizationGap > m_OptimizationGapThreshold
        && m_NumberCurrentSupportVectors < float( 0.9) * m_TrainingData->GetSize();
    }

    //! @brief finalize support vector model and determine final SVs, Alphas and Bias
    void ApproximatorSupportVectorMachine::FinalizeSupportVectorModel()
    {
      // reconstruct the whole gradient
      ReconstructGradient();

      // resort alpha vector to original length
      size_t alpha_count( 0);

      // counter for bound alphas
      size_t bound_alpha_count( 0);

      // counter for iteration
      size_t vector_counter( 0);

      // final vector with all alphas of support vectors
      storage::Vector< float> alpha_final;

      storage::List< size_t> sv_indices;

      // assembly of the final alpha vector of the lagrange multipliers
      for
      (
        storage::Vector< float>::const_iterator iter_alpha_begin( m_Alpha.Begin()),
        iter_alpha_middle( m_Alpha.Begin() + m_TrainingData->GetSize()),
        iter_alpha_middle_const( m_Alpha.Begin() + m_TrainingData->GetSize());
        iter_alpha_begin != iter_alpha_middle_const;
        ++iter_alpha_begin, ++iter_alpha_middle
      )
      {
        // individual lagrange multiplier value alpha
        // creating final alphas
        const float alpha_value( *iter_alpha_begin - *iter_alpha_middle);

        // if alpha is not zero then correspondent vector is a support vector
        if( alpha_value != 0)
        {
          alpha_final.PushBack( alpha_value);

          sv_indices.PushBack( vector_counter);

          ++alpha_count;
        }

        // increase iteration counter
        ++vector_counter;

        // if alpha lays has a value of penalty parameter C (boundary)
        if( alpha_value == m_CostParameterC)
        {
          ++bound_alpha_count;
        }
      }

      // matrix containing support vectors
      linal::Matrix< float> sv_matrix
      (
        sv_indices.GetSize(),                             // rows
        m_TrainingData->GetFeatureSize(),                 // cols
        float( 0)                                         // default value
      );

      // reset counter
      vector_counter = 0;

      // fill support vector matrix with vectors by support vector index
      for
      (
        storage::List< size_t>::const_iterator itr_sv( sv_indices.Begin()), itr_sv_end( sv_indices.End());
        itr_sv != itr_sv_end;
        ++itr_sv, ++vector_counter
      )
      {
        // get current support vector
        const FeatureReference< float> &support_vector( m_TrainingData->GetFeaturesPtr()->operator ()( *itr_sv));
        // copy support vector to support vector matrix
        std::copy( support_vector.Begin(), support_vector.End(), sv_matrix[ vector_counter]);
      }

      // assemble SV Model
      m_Model->SetAlpha( alpha_final);
      m_Model->SetBias( CalculateBias());
      m_Model->SetSupportVectors( FeatureDataSet< float>( sv_matrix));
      m_Model->SetNumberSupportVectors( m_Model->GetSupportVectors().GetNumberFeatures());
      m_Model->SetNumberBoundSupportVectors( bound_alpha_count);
    }

    //! read NeuralNetwork from std::istream
    std::istream &ApproximatorSupportVectorMachine::Read( std::istream &ISTREAM)
    {
      storage::Vector< float>::s_Instance.IsDefined();

      // read members
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
      io::Serialize::Read( m_NumberIterations, ISTREAM);
      io::Serialize::Read( m_OptimizationGapThreshold, ISTREAM);
      io::Serialize::Read( m_OptimizationGap, ISTREAM);
      io::Serialize::Read( m_LastObjectiveFunctionValue, ISTREAM);

      // return
      return ISTREAM;
    }

    //! write NeuralNetwork into std::ostream
    std::ostream &ApproximatorSupportVectorMachine::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_CostParameterC, OSTREAM, INDENT) << std::endl;
      io::Serialize::Write( m_Status, OSTREAM, INDENT) << std::endl;
      io::Serialize::Write( m_Alpha, OSTREAM, INDENT) << std::endl;
      io::Serialize::Write( m_Gradient, OSTREAM, INDENT) << std::endl;
      io::Serialize::Write( m_GradientBar, OSTREAM, INDENT) << std::endl;
      io::Serialize::Write( m_Bias, OSTREAM, INDENT) << std::endl;
      io::Serialize::Write( m_Signs, OSTREAM, INDENT) << std::endl;
      io::Serialize::Write( m_Labels, OSTREAM, INDENT) << std::endl;
      io::Serialize::Write( m_ActiveSize, OSTREAM, INDENT) << std::endl;
      io::Serialize::Write( m_ProbLength, OSTREAM, INDENT) << std::endl;
      io::Serialize::Write( m_Model, OSTREAM, INDENT) << std::endl;
      io::Serialize::Write( m_NumberIterations, OSTREAM, INDENT) << std::endl;
      io::Serialize::Write( m_OptimizationGapThreshold, OSTREAM, INDENT) << std::endl;
      io::Serialize::Write( m_OptimizationGap, OSTREAM, INDENT) << std::endl;
      io::Serialize::Write( m_LastObjectiveFunctionValue, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

    //! @brief method checks whether a LaGrange Multiplier reached a certain boundary or not
    //! @param ALPHA a LaGrange Multiplier
    //! @return const int that indicates whether ALPHA reached a certain boundary or not
    int ApproximatorSupportVectorMachine::AlphaToStatus( const float &ALPHA) const
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
    void ApproximatorSupportVectorMachine::UpdateAlphaStatus( const int &FEATURE_VECTOR_I)
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
    inline bool ApproximatorSupportVectorMachine::IsUpperBound( const size_t &FEATURE_VECTOR_I) const
    {
      return m_Status( FEATURE_VECTOR_I) == s_Upper_Bound;
    }

    //! @brief checks whether a feature_vector i reached the lower boundary
    //! @param FEATURE_VECTOR_I a position of a feature vector in m_Status
    //! @return const bool that indicates whether FEATURE_VECTOR_I reached lower boundary
    inline bool ApproximatorSupportVectorMachine::IsLowerBound( const size_t &FEATURE_VECTOR_I) const
    {
      return m_Status( FEATURE_VECTOR_I) == s_Lower_Bound;
    }

    //! @brief checks whether a feature_vector i is not bound and between the boundaries
    //! @param FEATURE_VECTOR_I a position of a feature vector in m_Status
    //! @return const bool that indicates whether FEATURE_VECTOR_I reached no boundary
    inline bool ApproximatorSupportVectorMachine::IsFree( const size_t &FEATURE_VECTOR_I) const
    {
      return m_Status( FEATURE_VECTOR_I) == s_Free;
    }

    //! @brief Initialize member vectors and variables to prepare SVR training
    void ApproximatorSupportVectorMachine::InitializeMemberVectorsForTraining()
    {
      // Initialization of variables and vectors
      const size_t number_feature_vectors( m_TrainingData->GetSize());

      // has to be set in case of this class is read from file and constructor was not explicitely applied
      m_ActiveSize = 0;
      m_ProbLength = 2 * number_feature_vectors;

      // initialize vector for lagrange multipliers
      m_Alpha = storage::Vector< float>( m_ProbLength, 0.0);

      // initialize vector for labeling of
      m_Labels = storage::Vector< float>();
      m_Labels.InsertElements( 0, 1.0, number_feature_vectors);
      m_Labels.InsertElements( number_feature_vectors, -1.0, number_feature_vectors);

      // initialize vector for signs of every label
      m_Signs = storage::Vector< int>();
      m_Signs.InsertElements( 0, 1, number_feature_vectors);
      m_Signs.InsertElements( number_feature_vectors, -1, number_feature_vectors);

      // status indication of KKT determination for every vector
      m_Status = storage::Vector< size_t>( m_ProbLength, size_t( 0));

      //
      m_GradientBar = storage::Vector< float>( m_ProbLength, 0.0);

      // initialization of bias for support vector model
      m_Bias = storage::Vector< float>( 2 * number_feature_vectors, m_P);

      // initialization of gradients for support vector model
      m_Gradient = storage::Vector< float>( 2 * number_feature_vectors);

      // initializing
      for( size_t progress( 0); progress < number_feature_vectors; ++progress)
      {
        const float result( m_TrainingData->GetResultsPtr()->operator ()( progress)( 0));
        m_Bias( progress) -= result;
        m_Bias( number_feature_vectors + progress) += result;
      }

      m_Gradient = m_Bias;

      // update alpha_status
      for( size_t progress( 0); progress < m_ProbLength; ++progress)
      {
        UpdateAlphaStatus( progress);
      }

      m_ActiveSize = m_ProbLength;

      linal::Vector< float> Q_i;

      // initializing Gradients
      for( size_t progress( 0); progress < m_ActiveSize; ++progress)
      {
        if( !IsLowerBound( progress))
        {
          Q_i = m_Model->GetKernel()->GetInputVectorIKernelMatrix( *( m_TrainingData->GetFeaturesPtr()), progress, m_Signs, m_ProbLength);

          for( size_t progress_internal( 0); progress_internal < m_ActiveSize; ++progress_internal)
          {
            m_Gradient( progress_internal) += m_Alpha( progress) * Q_i( progress_internal);
          }

          if( IsUpperBound( progress))
          {
            for( size_t progress_internal( 0); progress_internal < m_ActiveSize; ++progress_internal)
            {
              m_GradientBar( progress_internal) += m_CostParameterC * Q_i( progress_internal);
            }
          }
        }
      }

      BCL_MessageDbg( " Initialization.. complete");
    }

    //! @brief determine the bias value for Support Vector Regression Model
    //! @return calculated bias value
    float ApproximatorSupportVectorMachine::CalculateBias() const
    {
      float bias;
      int nr_free( 0);
      float sum_free( 0);

      // initialize upper and lower boundaries
      float upper_bound(  std::numeric_limits< float>::infinity());
      float lower_bound( -std::numeric_limits< float>::infinity());

      size_t progress( 0);
      // iterate over all labes and gradients to update upper and lower boundaries
      for
      (
        storage::Vector< float>::const_iterator iter_begin_labels( m_Labels.Begin()),
        iter_begin_gradient( m_Gradient.Begin());
        iter_begin_labels != m_Labels.End();
        ++iter_begin_labels, ++iter_begin_gradient, ++progress
      )
      {
        // compute gradient of particular label
        const float label_gradient( *iter_begin_labels * ( *iter_begin_gradient));

        // check for boundary conditions and update upper or lower boundaries accordingly
        if( IsLowerBound( progress))
        {
          if( m_Labels( progress) > 0)
          {
            upper_bound = std::min( upper_bound, label_gradient);
          }
          else
          {
            lower_bound = std::max( lower_bound, label_gradient);
          }
        }
        else if( IsUpperBound( progress))
        {
          if( m_Labels( progress) < 0)
          {
            upper_bound = std::min( upper_bound, label_gradient);
          }
          else
          {
            lower_bound = std::max( lower_bound, label_gradient);
          }
        }
        else
        {
          // case where no boundary was hit
          ++nr_free;
          sum_free += label_gradient;
        }
      }

      // calculate bias
      if( nr_free > 0)
      {
        bias = sum_free / nr_free;
      }
      else
      {
        bias = ( upper_bound + lower_bound) / 2;
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
    float ApproximatorSupportVectorMachine::DetermineFeatureVectorCombination
    (
      size_t &FIRST_VECTOR_INDEX,
      size_t &SECOND_VECTOR_INDEX
    ) const
    {
      // return i,j such that
      // i: maximizes -y_i * grad(f)_i, i in I_up(\alpha)
      // j: minimizes the decrease of obj value
      //    (if quadratic coefficient <= 0, replace it with tau)
      //    -y_j*grad(f)_j < -y_i*grad(f)_i, j in I_low(\alpha)

      float maximum_gradient( -1 * std::numeric_limits< float>::infinity()); //-INFINITY;
      float Gmax2( -1 * std::numeric_limits< float>::infinity()); //-INFINITY;

      int maximum_gradient_index( -1);
      int Gmin_idx( -1);

      size_t process( 0);

      for
      (
        storage::Vector< float>::const_iterator
          iter_begin_gradient( m_Gradient.Begin()),
          iter_end_gradient( m_Gradient.End()),
          iter_begin_labels( m_Labels.Begin());
        iter_begin_gradient != iter_end_gradient;
        ++iter_begin_gradient, ++iter_begin_labels, ++process
      )
      {
        if( *iter_begin_labels == +1)
        {
          if( !IsUpperBound( process))
          {
            if( -( *iter_begin_gradient) >= maximum_gradient)
            {
              maximum_gradient = -( *iter_begin_gradient);
              maximum_gradient_index = process;
            }
          }
        }
        else
        {
          if( !IsLowerBound( process))
          {
            if( *iter_begin_gradient >= maximum_gradient)
            {
              maximum_gradient = *iter_begin_gradient;
              maximum_gradient_index = process;
            }
          }
        }
      }

      const int index_i( maximum_gradient_index);

      linal::Vector< float> Q_i;

      if( index_i != -1) // null Q_i not accessed: Gmax=-INF if index_i=-1
      {
        Q_i = m_Model->GetKernel()->GetInputVectorIKernelMatrix
        (
          *( m_TrainingData->GetFeaturesPtr()),
          index_i,
          m_Signs,
          m_ProbLength
        );
      }

      linal::Vector< float> QD( m_ProbLength, float( 0.0));

      float kernel_value( 0.0);
      size_t progress( 0);

      for
      (
        linal::Vector< float>::iterator iter_begin_QD( QD.Begin()),
        iter_middle_QD( QD.Begin() + QD.GetSize() / 2),
        iter_end_QD( QD.End());
        iter_middle_QD != iter_end_QD;
        ++iter_begin_QD, ++iter_middle_QD, ++progress
      )
      {

        const FeatureReference< float> &item_ref_at_position( m_TrainingData->GetFeaturesPtr()->operator ()( progress));

        kernel_value = m_Model->GetKernel()->operator()
        (
          item_ref_at_position,
          item_ref_at_position
        );

        *iter_begin_QD  = kernel_value;
        *iter_middle_QD = kernel_value;
      }

      float obj_diff_min( std::numeric_limits< float>::infinity()); //+INFINITY;

      float grad_diff( 0);
      float obj_diff( 0);
      float quad_coef( 0);

      size_t index_j( 0);

      linal::Vector< float>::const_iterator iter_begin_q_i( Q_i.Begin());

      for
      (
        storage::Vector< float>::const_iterator
          iter_begin_labels( m_Labels.Begin()),
          iter_end_labels( m_Labels.End()),
          iter_begin_gradient( m_Gradient.Begin());
        iter_begin_labels != iter_end_labels;
        ++iter_begin_labels, ++iter_begin_gradient, ++iter_begin_q_i, ++index_j
      )
      {
        // if label is positive
        if( *iter_begin_labels == +1)
        {
          if( !IsLowerBound( index_j))
          {
            grad_diff = maximum_gradient + *iter_begin_gradient;

            if( -*iter_begin_gradient >= Gmax2)
            {
              Gmax2 = *iter_begin_gradient;
            }

            if( grad_diff > 0)
            {
              // QD (k) = QD ( k + l)
              quad_coef = Q_i( index_i) + QD( index_j) - 2 * m_Labels( index_i) * ( *iter_begin_q_i);

              if( quad_coef > 0)
              {
                obj_diff = -( grad_diff * grad_diff) / quad_coef;
              }
              else
              {
                obj_diff = -( grad_diff * grad_diff) / m_EPS_A;
              }

              if( obj_diff <= obj_diff_min)
              {
                Gmin_idx = index_j;
                obj_diff_min = obj_diff;
              }
            }
          }
        }
        else
        {
          // if label is negative
          if( !IsUpperBound( index_j))
          {
            grad_diff = maximum_gradient - *iter_begin_gradient;

            if( -*iter_begin_gradient >= Gmax2)
            {
              Gmax2 = -( *iter_begin_gradient);
            }

            if( grad_diff > 0)
            {
              quad_coef = Q_i( index_i) + QD( index_j) + 2 * m_Labels( index_i) * ( *iter_begin_q_i);

              if( quad_coef > 0)
              {
                obj_diff = -( grad_diff * grad_diff) / quad_coef;
              }
              else
              {
                obj_diff = -( grad_diff * grad_diff) / m_EPS_A;
              }

              if( obj_diff <= obj_diff_min)
              {
                Gmin_idx = index_j;
                obj_diff_min = obj_diff;
              }
            }
          }
        }
      }

      FIRST_VECTOR_INDEX  = maximum_gradient_index;
      SECOND_VECTOR_INDEX = Gmin_idx;

      return maximum_gradient + Gmax2; // optimization gap
    }

    //! @brief solve sub problem of quadratic problem according two feature vectors
    //! @param SVR_MODEL support vector machine model of interest
    //! @param FEATURE_VECTOR_I index of feature vector in MATRIX
    //! @param FEATURE_VECTOR_J index of feature vector in MATRIX
    //! @return a pair of vector of float with computed values of kernel function
    void ApproximatorSupportVectorMachine::SolveQuadraticProblemSubProblem
    (
      const int &FEATURE_VECTOR_I,
      const int &FEATURE_VECTOR_J
    )
    {
      const float old_m_alpha_i( m_Alpha( FEATURE_VECTOR_I));
      const float old_m_alpha_j( m_Alpha( FEATURE_VECTOR_J));

      // using Kernel Q_i(x) = K(i,x)
      linal::Vector< float> Q_i
      (
        m_Model->GetKernel()->GetInputVectorIKernelMatrix( *( m_TrainingData->GetFeaturesPtr()), FEATURE_VECTOR_I, m_Signs, m_ProbLength)
      );

      // using Kernel Q_j(x) = K(j,x)
      linal::Vector< float> Q_j
      (
        m_Model->GetKernel()->GetInputVectorIKernelMatrix( *( m_TrainingData->GetFeaturesPtr()), FEATURE_VECTOR_J, m_Signs, m_ProbLength)
      );

      float quad_coef( Q_i( FEATURE_VECTOR_I) + Q_j( FEATURE_VECTOR_J) + 2 * Q_i( FEATURE_VECTOR_J));

      if( quad_coef <= 0)
      {
        quad_coef = m_EPS_A;
      }

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

      auto iter_begin_Qi( Q_i.Begin()), iter_begin_Qj( Q_j.Begin());
      for
      (
        storage::Vector< float>::iterator
          iter_begin_gradient( m_Gradient.Begin()),
          iter_end_gradient( m_Gradient.End());
        iter_begin_gradient != iter_end_gradient;
        ++iter_begin_gradient, ++iter_begin_Qi, ++iter_begin_Qj
      )
      {
        *iter_begin_gradient += *iter_begin_Qi * delta_alpha_i + *iter_begin_Qj * delta_alpha_j;
      }

      // update alpha_status and m_GradientBar
      const bool feature_i_is_upper_bound( IsUpperBound( FEATURE_VECTOR_I));
      const bool feature_j_is_upper_bound( IsUpperBound( FEATURE_VECTOR_J));

      UpdateAlphaStatus( FEATURE_VECTOR_I);
      UpdateAlphaStatus( FEATURE_VECTOR_J);

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
        iter_begin_Qi = Q_i.Begin();
        const double effective_cost( cost_multiplier * m_CostParameterC);
        for
        (
          storage::Vector< float>::iterator
            iter_begin_gradient_bar( m_GradientBar.Begin()),
            iter_end_gradient_bar( m_GradientBar.End());
          iter_begin_gradient_bar != iter_end_gradient_bar;
          ++iter_begin_gradient_bar, ++iter_begin_Qi
        )
        {
          *iter_begin_gradient_bar += effective_cost * ( *iter_begin_Qi);
        }
      }

    }

    //! @brief reconstruct the gradients for the final model
    void ApproximatorSupportVectorMachine::ReconstructGradient()
    {
      // reconstruct inactive elements of m_Gradient from m_GradientBar and free variables
      if( m_ActiveSize == m_ProbLength)
      {
        return;
      }

      for
      (
        storage::Vector< float>::iterator iter_begin_gradient( m_Gradient.Begin()),
        iter_end_gradient( m_Gradient.End()),
        iter_begin_gradient_bar( m_GradientBar.Begin()),
        iter_begin_bias( m_Bias.Begin());
        iter_begin_gradient != iter_end_gradient;
        ++iter_begin_gradient, ++iter_begin_gradient_bar, ++iter_begin_bias
      )
      {
        *iter_begin_gradient = *iter_begin_gradient_bar + *iter_begin_bias;
      }

      // vector of kernel values of vector i and all other vectors in training set
      linal::Vector< float> Q_i;

      // progress counter accessing elements at a certain position in vector
      size_t progress( 0);

      // iterate over all LaGrange multipliers
      for
      (
        storage::Vector< float>::iterator
          iter_begin_alpha( m_Alpha.Begin()),
          iter_end_alpha( m_Alpha.End());
        iter_begin_alpha != iter_end_alpha;
        ++iter_begin_alpha, ++progress
      )
      {
        // if correspondent status values indicates that it is not between upper or lower bound
        if( IsFree( progress))
        {
          // compute kernel vector for specific input vector i
          Q_i = m_Model->GetKernel()->GetInputVectorIKernelMatrix( *( m_TrainingData->GetFeaturesPtr()), progress, m_Signs, m_ProbLength);

          linal::Vector< float>::iterator iter_begin_Qi( Q_i.Begin());

          // iterate over all feature vector gradients
          for
          (
            storage::Vector< float>::iterator
              iter_begin_gradient( m_Gradient.Begin() + m_ActiveSize),
              iter_end_gradient( m_Gradient.End());
            iter_begin_gradient != iter_end_gradient;
            ++iter_begin_gradient, ++iter_begin_Qi
          )
          {
            // add next alpha * kernel(i,i) to gradient
            *iter_begin_gradient += *iter_begin_alpha * ( *iter_begin_Qi);
          }
        }
      }
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer ApproximatorSupportVectorMachine::GetSerializer() const
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
      return parameters;
    }

  } // namespace model
} // namespace bcl
