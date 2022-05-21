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
#include "model/bcl_model_approximator_linear_regression.h"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_matrix_inversion_cholesky.h"
#include "linal/bcl_linal_matrix_inversion_moore_penrose.h"
#include "util/bcl_util_stopwatch.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {
    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> ApproximatorLinearRegression::s_Instance
    (
      util::Enumerated< ApproximatorBase>::AddInstance( new ApproximatorLinearRegression())
    );

    //! @brief default constructor
    ApproximatorLinearRegression::ApproximatorLinearRegression() :
      m_MatrixSolver( linal::MatrixInversionCholesky< float>()),
      m_LastObjectiveFunctionResult( util::GetUndefined< float>())
    {
    }

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ApproximatorLinearRegression::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &ApproximatorLinearRegression::GetAlias() const
    {
      static const std::string s_alias( "LinearRegression");
      return s_alias;
    }

    //! @brief set training data set for a specific iterate in approximator framework
    //! @param DATA training data set
    void ApproximatorLinearRegression::SetTrainingData
    (
      util::ShPtr< descriptor::Dataset> &DATA
    )
    {
      BCL_Assert( !DATA->IsEmpty(), "no training data loaded");

      // rescale function for in an output and denormalization
      m_TrainingData = DATA;
      linal::Matrix< float> coefficients( m_TrainingData->GetFeatureSize(), m_TrainingData->GetResultSize());
      m_Weights = util::ShPtr< Interface>( new MultipleLinearRegression( coefficients));
      BCL_MessageStd( "Setting up training data with " + util::Format()( DATA->GetSize()) + " points");
    }

    //! @brief construct a model from the current iterate
    //! @return shptr to the new model interface
    util::ShPtr< Interface> ApproximatorLinearRegression::GetCurrentModel() const
    {
      return m_Weights;
    }

    //! @brief returns the current approximation
    //! @return current argument result pair
    const util::ShPtr< storage::Pair< util::ShPtr< Interface>, float> >
      ApproximatorLinearRegression::GetCurrentApproximation() const
    {
      return
        util::ShPtr< storage::Pair< util::ShPtr< Interface>, float> >
        (
          new storage::Pair< util::ShPtr< Interface>, float>( m_Weights, m_LastObjectiveFunctionResult)
        );
    }

    //! @brief conducts the next approximation step and stores the approximation
    void ApproximatorLinearRegression::Next()
    {
      // check whether the model has already been trained
      if( !util::IsNaN( m_LastObjectiveFunctionResult))
      {
        return;
      }

      // assert when there is no interval
      BCL_Assert( m_TrainingData.IsDefined() && !m_TrainingData->IsEmpty(), "no training data given!");

      const size_t number_features( m_TrainingData->GetSize());
      const size_t feature_size( m_TrainingData->GetFeatureSize());
      const size_t result_size( m_TrainingData->GetResultSize());
      BCL_Assert( number_features >= feature_size, "insufficient training data to perform regression!");

      util::Implementation< linal::MatrixInversionInterface< float> > solver( m_MatrixSolver);

      const linal::MatrixConstReference< float> features( m_TrainingData->GetFeatures().GetMatrix());
      const linal::MatrixConstReference< float> results( m_TrainingData->GetResults().GetMatrix());

      // compute the product matrix
      linal::Matrix< float> transpose_matrix_times_matrix;
      linal::Matrix< float> transpose_matrix_times_results;
      linal::MatrixConstReference< float> matrix_to_solve( features);
      linal::MatrixConstReference< float> b_vectors( results);
      bool is_moore_penrose( m_MatrixSolver->GetAlias() == linal::MatrixInversionMoorePenrose< float>().GetAlias());
      // Moore penrose computes the matrix inverse directly, but for all other methods,
      // it is necessary to compute X^T * X (lhs) and X^T * y (rhs)
      if( number_features > feature_size && !is_moore_penrose)
      {
        BCL_MessageStd
        (
          "Multiplying matrices of size " + util::Format()( number_features)
          + " X " + util::Format()( feature_size)
        );
        util::Stopwatch mult_timer
        (
          "feature Transpose * features & feature Transpose * results",
          util::Time(),
          util::Message::e_Standard
        );
        transpose_matrix_times_matrix = linal::MatrixTransposeTimesMatrix( features);
        matrix_to_solve.Reference( transpose_matrix_times_matrix);
        transpose_matrix_times_results = linal::MatrixTransposeTimesMatrix( features, results);
        b_vectors.Reference( transpose_matrix_times_results); // feature_size X result_size
      }

      BCL_MessageStd( "Solving matrix of size " + util::Format()( matrix_to_solve.GetNumberRows()));

      util::Stopwatch set_matrix_timer( "Decomposition/Inversion", util::Time(), util::Message::e_Standard, false);
      // compute the inverse of the features matrix
      if( !solver->SetMatrix( matrix_to_solve))
      {
        set_matrix_timer.Stop();
        // let the user know that their chosen method of matrix solution failed, use moore-penrose as backup
        BCL_MessageStd
        (
          "Matrix solver " + solver->GetAlias() + " could not handle the matrix; defaulting to moore-penrose algorithm"
        );
        solver = linal::MatrixInversionMoorePenrose< float>();
        // reset the matrix times matrix transpose to conserve memory
        transpose_matrix_times_matrix = linal::Matrix< float>();
        solver->SetMatrix( features);
        is_moore_penrose = true;
      }
      set_matrix_timer.Stop();
      set_matrix_timer.WriteMessage();

      util::Stopwatch solve_matrix_timer( "Solve Matrix", util::Time(), util::Message::e_Standard);
      linal::Matrix< float> coefficients;
      if( is_moore_penrose)
      {
        // compute the coefficients
        coefficients = solver->ComputeInverse() * results;
      }
      else
      {
        // solve for the coefficients
        coefficients = b_vectors;
        for( size_t result_num( 0); result_num < result_size; ++result_num)
        {
          coefficients.ReplaceCol( result_num, solver->Solve( coefficients.GetCol( result_num)));
        }
      }
      solve_matrix_timer.Stop();

      BCL_MessageStd( "Calculated coefficients");

      // create the model
      m_Weights = util::ShPtr< Interface>( new MultipleLinearRegression( coefficients));

      // store the objective function result for later iterations
      m_LastObjectiveFunctionResult = m_ObjectiveFunction->operator ()( m_Weights);
      BCL_MessageStd
      (
        " ObjFunction: " + util::Format()( m_LastObjectiveFunctionResult)
      );

      this->GetTracker().Track( GetCurrentApproximation());
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ApproximatorLinearRegression::Read( std::istream &ISTREAM)
    {
      ApproximatorBase::Read( ISTREAM);
      io::Serialize::Read( m_Weights, ISTREAM);
      io::Serialize::Read( m_LastObjectiveFunctionResult, ISTREAM);
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &ApproximatorLinearRegression::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      ApproximatorBase::Write( OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Weights, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_LastObjectiveFunctionResult, OSTREAM, INDENT);
      return OSTREAM;
    }

    //! @brief evaluates whether the approximation can continue
    //! @return true, if the approximation can continue - otherwise false
    bool ApproximatorLinearRegression::CanContinue() const
    {
      // this approximation cannot be improved as it is analytic
      return !this->GetTracker().GetIteration();
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer ApproximatorLinearRegression::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Performs multiple linear regression see http://en.wikipedia.org/wiki/Linear_regression"
      );

      parameters.AddInitializer
      (
        "objective function",
        "function that evaluates the model after each iteration",
        io::Serialization::GetAgent( &m_ObjectiveFunction->GetImplementation()),
        "RMSD"
      );
      parameters.AddInitializer
      (
        "solver",
        "Solution method for the matrix",
        io::Serialization::GetAgent( &m_MatrixSolver),
        "Cholesky"
      );
      return parameters;
    }

  } // namespace model
} // namespace bcl
