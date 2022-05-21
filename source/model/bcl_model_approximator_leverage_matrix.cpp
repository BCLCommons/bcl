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
#include "model/bcl_model_approximator_leverage_matrix.h"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_matrix_inversion_cholesky.h"
#include "linal/bcl_linal_matrix_inversion_moore_penrose.h"
#include "model/bcl_model_leverage_matrix.h"
#include "util/bcl_util_stopwatch.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {
    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> ApproximatorLeverageMatrix::s_Instance
    (
      util::Enumerated< ApproximatorBase>::AddInstance( new ApproximatorLeverageMatrix())
    );

    //! @brief default constructor
    ApproximatorLeverageMatrix::ApproximatorLeverageMatrix() :
      m_LastObjectiveFunctionResult( util::GetUndefined< float>())
    {
    }

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ApproximatorLeverageMatrix::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &ApproximatorLeverageMatrix::GetAlias() const
    {
      static const std::string s_alias( "Leverage");
      return s_alias;
    }

    //! @brief set training data set for a specific iterate in approximator framework
    //! @param DATA training data set
    void ApproximatorLeverageMatrix::SetTrainingData
    (
      util::ShPtr< descriptor::Dataset> &DATA
    )
    {
      BCL_Assert( !DATA->IsEmpty(), "no training data loaded");

      // rescale function for in an output and denormalization
      m_TrainingData = DATA;
      m_TrainingData->GetFeatures().Rescale( math::Range< float>( 0.0, 1.0), RescaleFeatureDataSet::e_AveStd);

      linal::Matrix< float> coefficients( m_TrainingData->GetFeatureSize(), m_TrainingData->GetFeatureSize());
      double av_hat( double( m_TrainingData->GetFeatureSize() + 1) / double( m_TrainingData->GetSize()));
      coefficients /= float( av_hat);
      m_LeverageMatrix = util::ShPtr< Interface>( new LeverageMatrix( coefficients, GetRescaleFeatureDataSet()));
      BCL_MessageStd( "Setting up training data with " + util::Format()( DATA->GetSize()) + " points");
    }

    //! @brief construct a model from the current iterate
    //! @return shptr to the new model interface
    util::ShPtr< Interface> ApproximatorLeverageMatrix::GetCurrentModel() const
    {
      return m_LeverageMatrix;
    }

    //! @brief returns the current approximation
    //! @return current argument result pair
    const util::ShPtr< storage::Pair< util::ShPtr< Interface>, float> >
      ApproximatorLeverageMatrix::GetCurrentApproximation() const
    {
      return
        util::ShPtr< storage::Pair< util::ShPtr< Interface>, float> >
        (
          new storage::Pair< util::ShPtr< Interface>, float>( m_LeverageMatrix, m_LastObjectiveFunctionResult)
        );
    }

    //! @brief conducts the next approximation step and stores the approximation
    void ApproximatorLeverageMatrix::Next()
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
      BCL_Assert( number_features >= feature_size, "insufficient training data to perform regression!");

      const linal::MatrixConstReference< float> features( m_TrainingData->GetFeatures().GetMatrix());
      const linal::MatrixConstReference< float> results( m_TrainingData->GetResults().GetMatrix());

      // compute the product matrix
      linal::Matrix< float> transpose_matrix_times_matrix( linal::MatrixTransposeTimesMatrix( features));

      BCL_MessageStd( "Solving matrix of size " + util::Format()( transpose_matrix_times_matrix.GetNumberRows()));
      linal::MatrixInversionMoorePenrose< float> mp_solver;
      // compute the inverse of the features matrix
      mp_solver.SetMatrix( transpose_matrix_times_matrix);

      const float av_hat( double( m_TrainingData->GetFeatureSize() + 1) / double( m_TrainingData->GetSize()));
      linal::Matrix< float> mp_inverse( mp_solver.ComputeInverse());
      mp_inverse /= av_hat;
      // create the model
      m_LeverageMatrix = util::ShPtr< Interface>( new LeverageMatrix( mp_inverse, GetRescaleFeatureDataSet()));

      // store the objective function result for later iterations
      m_LastObjectiveFunctionResult = 0.0;

      this->GetTracker().Track( GetCurrentApproximation());
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ApproximatorLeverageMatrix::Read( std::istream &ISTREAM)
    {
      ApproximatorBase::Read( ISTREAM);
      io::Serialize::Read( m_LeverageMatrix, ISTREAM);
      io::Serialize::Read( m_LastObjectiveFunctionResult, ISTREAM);
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &ApproximatorLeverageMatrix::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      ApproximatorBase::Write( OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_LeverageMatrix, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_LastObjectiveFunctionResult, OSTREAM, INDENT);
      return OSTREAM;
    }

    //! @brief evaluates whether the approximation can continue
    //! @return true, if the approximation can continue - otherwise false
    bool ApproximatorLeverageMatrix::CanContinue() const
    {
      // this approximation cannot be improved as it is analytic
      return !this->GetTracker().GetIteration();
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer ApproximatorLeverageMatrix::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Computes the leverage matrix (aka projection or hat matrix), which allows identification of significant outliers "
        "that would likely substantially influence any simple linear model of the system. A returned value > 2 represents "
        "probable outliers, while greater than 3 represent definitive outliers. The average value is 1 for all values in "
        "the training set"
      );

      parameters.AddInitializer
      (
        "objective function",
        "function that evaluates the model after each iteration",
        io::Serialization::GetAgent( &m_ObjectiveFunction->GetImplementation()),
        "RMSD"
      );
      return parameters;
    }

  } // namespace model
} // namespace bcl
