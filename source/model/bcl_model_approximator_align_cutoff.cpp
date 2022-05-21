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
#include "model/bcl_model_approximator_align_cutoff.h"

// includes from bcl - sorted by alphabetically
#include "math/bcl_math_contingency_matrix.h"
#include "math/bcl_math_roc_curve.h"
#include "model/bcl_model_align_cutoff.h"
#include "model/bcl_model_objective_function_constant.h"

// external includes - sorted by alphabetically

namespace bcl
{
  namespace model
  {

    // add the interface to the set of known implementations
    const util::SiPtr< const util::ObjectInterface> ApproximatorAlignCutoff::s_Instance
    (
      util::Enumerated< ApproximatorBase>::AddInstance( new ApproximatorAlignCutoff())
    );

    //! @brief get the numerical tolerance used in this class
    //! @return the numerical tolerance used in this class
    const float &ApproximatorAlignCutoff::GetNumericalTolerance()
    {
      static const float s_tolerance( 1.0e-5);
      return s_tolerance;
    }

    //! @brief default constructor
    ApproximatorAlignCutoff::ApproximatorAlignCutoff() :
      m_InternalIterate(),
      m_DescaledExperimentalCutoff( 0.5),
      m_ExperimentalCutoff( 0.5),
      m_ModelCutoff( 0.5),
      m_ExperimentalAverageAboveCutoff( 1.0),
      m_ExperimentalAverageBelowCutoff( 0.0),
      m_ModelAverageAboveCutoff( 1.0),
      m_ModelAverageBelowCutoff( 0.0),
      m_PositivesAboveThreshold( false),
      m_NumberActives( 0)
    {
    }

    //! @brief copy constructor
    //! @return a new ApproximatorAlignCutoff copied from this instance
    ApproximatorAlignCutoff *ApproximatorAlignCutoff::Clone() const
    {
      return new ApproximatorAlignCutoff( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ApproximatorAlignCutoff::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &ApproximatorAlignCutoff::GetAlias() const
    {
      static const std::string s_Name( "AlignCutoff");
      return s_Name;
    }

    //! @brief construct a model from the current iterate
    //! @return shptr to the new model interface
    util::ShPtr< Interface> ApproximatorAlignCutoff::GetCurrentModel() const
    {
      // make a new model out of the current data members
      return
        util::ShPtr< Interface>
        (
          new AlignCutoff
          (
            m_LastInternallyGeneratedModel.IsDefined() ? m_LastInternallyGeneratedModel : m_InternalIterate->GetCurrentModel(),
            GetRescaleFeatureDataSet(),
            GetRescaleResultDataSet(),
            m_ExperimentalCutoff,
            m_ModelCutoff,
            ( m_ExperimentalAverageAboveCutoff - m_ExperimentalCutoff) / ( m_ModelAverageAboveCutoff - m_ModelCutoff),
            ( m_ExperimentalAverageBelowCutoff - m_ExperimentalCutoff) / ( m_ModelAverageBelowCutoff - m_ModelCutoff)
          )
        );
    }

    //! @brief returns the current approximation
    //! @return current argument result pair
    const util::ShPtr< storage::Pair< util::ShPtr< Interface>, float> >
      ApproximatorAlignCutoff::GetCurrentApproximation() const
    {
      return
        util::ShPtr< storage::Pair< util::ShPtr< Interface>, float> >
        (
          new storage::Pair< util::ShPtr< Interface>, float>
          (
            m_LastModel,
            m_LastObjectiveFunctionValue
          )
        );
    }

    //! @brief set objective function to evaluate a monitoring dataset
    //! @param OBJ objective function of interest
    void ApproximatorAlignCutoff::SetObjectiveFunction( const util::ShPtr< ObjectiveFunctionWrapper> &OBJ)
    {
      ApproximatorBase::SetObjectiveFunction( OBJ);

      // reset the last internally generated model to force re-evaluation of the objective function if it is
      // called again
      m_LastInternallyGeneratedModel = util::ShPtr< Interface>();
    }

    //! @brief set training data set for a specific iterate in approximater framework
    //! @param DATA training data set
    void ApproximatorAlignCutoff::SetTrainingData
    (
      util::ShPtr< descriptor::Dataset> &DATA
    )
    {
      m_TrainingData = DATA;
      m_InternalIterate->SetTrainingData( DATA);
      const FeatureDataSetInterface< float> &results( *m_InternalIterate->GetTrainingData()->GetResultsPtr());

      m_ExperimentalCutoff = m_DescaledExperimentalCutoff;
      if( m_InternalIterate->GetRescaleResultDataSet().IsDefined())
      {
        // need to adjust the experimental threshold
        m_ExperimentalCutoff = m_InternalIterate->GetRescaleResultDataSet()->RescaleValue( 0, m_DescaledExperimentalCutoff);
      }

      math::RunningAverage< float> average_above_cutoff, average_below_cutoff;
      for( size_t row_id( 0), num_rows( results.GetNumberFeatures()); row_id < num_rows; ++row_id)
      {
        const float &result_value( results[ row_id][ 0]);

        // average values on each side of the cutoff cutoff
        if( result_value > m_ExperimentalCutoff)
        {
          average_above_cutoff += result_value;
        }
        else
        {
          average_below_cutoff += result_value;
        }
      }
      m_ExperimentalAverageAboveCutoff = average_above_cutoff.GetWeight() ? average_above_cutoff : m_ExperimentalCutoff;
      m_ExperimentalAverageBelowCutoff = average_below_cutoff.GetWeight() ? average_below_cutoff : m_ExperimentalCutoff;
      if( math::EqualWithinTolerance( m_ExperimentalAverageAboveCutoff, m_ExperimentalCutoff, GetNumericalTolerance()))
      {
        // essentially everything above the cutoff shares the same value
        // To allow ranking, and to avoid divisions by zero, make the average above cutoff slightly larger
        m_ExperimentalAverageAboveCutoff =
            m_ExperimentalCutoff > GetNumericalTolerance()
            ? m_ExperimentalCutoff * ( 1.0 + GetNumericalTolerance())
            : m_ExperimentalCutoff < -GetNumericalTolerance()
              ? m_ExperimentalCutoff / ( 1.0 + GetNumericalTolerance())
              : m_ExperimentalCutoff + GetNumericalTolerance();
      }

      if( math::EqualWithinTolerance( m_ExperimentalAverageBelowCutoff, m_ExperimentalCutoff, GetNumericalTolerance()))
      {
        // essentially everything below the cutoff shares the same value
        // To allow ranking, and to avoid divisions by zero, make the average above cutoff slightly larger
        m_ExperimentalAverageBelowCutoff =
            m_ExperimentalCutoff > GetNumericalTolerance()
            ? m_ExperimentalCutoff / ( 1.0 + GetNumericalTolerance())
            : m_ExperimentalCutoff < -GetNumericalTolerance()
              ? m_ExperimentalCutoff * ( 1.0 + GetNumericalTolerance())
              : m_ExperimentalCutoff - GetNumericalTolerance();
      }

      m_ModelAverageAboveCutoff = m_ExperimentalAverageAboveCutoff;
      m_ModelAverageBelowCutoff = m_ExperimentalAverageBelowCutoff;
      m_NumberActives = m_PositivesAboveThreshold ? average_above_cutoff.GetWeight() : average_below_cutoff.GetWeight();

      BCL_MessageVrb( "Initialize num actives: " + util::Format()( m_NumberActives));
      BCL_MessageVrb
      (
        "Initialize average above / below: " + util::Format()( m_ExperimentalAverageAboveCutoff)
        + " / " + util::Format()( m_ExperimentalAverageBelowCutoff)
      );

    } // SetTrainingData

    //! @brief conducts the next approximation step and stores the approximation
    void ApproximatorAlignCutoff::Next()
    {
      m_InternalIterate->Next();
      util::ShPtr< Interface> internal_model( m_InternalIterate->GetCurrentModel());

      // check whether the model was actually changed
      if( internal_model != m_LastInternallyGeneratedModel)
      {
        // model changed, recompute the cutoffs and create the new model
        m_LastInternallyGeneratedModel = internal_model;

        // determine the new cutoffs and averages
        DetermineAdjustedCutoff();

        // create a new classifier
        m_LastModel = GetCurrentModel();

        // compute the objective function
        m_LastObjectiveFunctionValue = ( *m_ObjectiveFunction)( m_LastModel);
      }

      // combine it with the objective function evaluation, best precision seen with associated cutoff
      this->GetTracker().Track( GetCurrentApproximation());
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer ApproximatorAlignCutoff::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Adjusts the values such that the model's output can be used with the input cutoff for rank-classification purposes "
        "Results of different machine learning methods or cross validation runs should not normally be averaged, for the "
        "purposes of making a jury prediction or presenting the results of a ranked-classification model, without using "
        "this wrapper.   This iterate applies a linear function to map the values given by the internal iterate.  The "
        "function is created so as to maximize precision of predictions, while maintaining a similar hit rate to the training data"
      );

      parameters.AddInitializer
      (
        "iterate",
        "iterate string that defines the machine learning training algorithm",
        io::Serialization::GetAgent( &m_InternalIterate)
      );

      parameters.AddInitializer
      (
        "cutoff",
        "potency cutoff",
        io::Serialization::GetAgent( &m_DescaledExperimentalCutoff)
      );

      parameters.AddInitializer
      (
        "parity",
        "specifies whether to prefer models the predict values above the cutoff (1) or below the cutoff (0)",
        io::Serialization::GetAgent( &m_PositivesAboveThreshold),
        "0"
      );

      return parameters;
    }

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERROR_STREAM the stream to write errors to
    bool ApproximatorAlignCutoff::ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM)
    {
      // acquire the objective function given to the internal iterate
      SetObjectiveFunction( m_InternalIterate->GetObjectiveFunction());

      // give the internal iterate an undefined objective function since its output will not be used
      m_InternalIterate->SetObjectiveFunction
      (
        util::ShPtr< ObjectiveFunctionWrapper>
        (
          new ObjectiveFunctionWrapper
          (
            ObjectiveFunctionConstant
            (
              0.0,
              m_ObjectiveFunction->GetImprovementType(),
              m_ObjectiveFunction->GetGoalType(),
              m_ObjectiveFunction->GetThreshold()
            )
          )
        )
      );

      // set the model cutoff == the experimental cutoff initially
      m_ExperimentalCutoff = m_DescaledExperimentalCutoff;
      m_ModelCutoff = m_ExperimentalCutoff;

      // set the experimental average below cutoff slightly below the threshold, to prevent divisions by zero
      m_ExperimentalAverageBelowCutoff =
        m_ExperimentalCutoff > 1.0
        ? m_ExperimentalCutoff / 2.0
        : m_ExperimentalCutoff < -1.0
          ? m_ExperimentalCutoff * 2.0
          : m_ExperimentalCutoff - 1.0;

      // set the experimental average above cutoff above the threshold, to prevent divisions by zero
      m_ExperimentalAverageAboveCutoff =
        m_ExperimentalCutoff > 1.0
        ? m_ExperimentalCutoff * 2.0
        : m_ExperimentalCutoff < -1.0
          ? m_ExperimentalCutoff / 2.0
          : m_ExperimentalCutoff + 1.0;

      m_ModelAverageAboveCutoff = m_ExperimentalAverageAboveCutoff;
      m_ModelAverageBelowCutoff = m_ExperimentalAverageBelowCutoff;

      return true;
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ApproximatorAlignCutoff::Read( std::istream &ISTREAM)
    {
      descriptor::Dataset::s_Instance.IsDefined();

      // read members
      io::Serialize::Read( m_InternalIterate, ISTREAM);
      io::Serialize::Read( m_DescaledExperimentalCutoff, ISTREAM);
      io::Serialize::Read( m_ExperimentalCutoff, ISTREAM);
      io::Serialize::Read( m_ModelCutoff, ISTREAM);
      io::Serialize::Read( m_ExperimentalAverageAboveCutoff, ISTREAM);
      io::Serialize::Read( m_ExperimentalAverageBelowCutoff, ISTREAM);
      io::Serialize::Read( m_ModelAverageAboveCutoff, ISTREAM);
      io::Serialize::Read( m_ModelAverageBelowCutoff, ISTREAM);
      io::Serialize::Read( m_PositivesAboveThreshold, ISTREAM);
      io::Serialize::Read( m_NumberActives, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &ApproximatorAlignCutoff::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_InternalIterate, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_DescaledExperimentalCutoff, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_ExperimentalCutoff, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_ModelCutoff, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_ExperimentalAverageAboveCutoff, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_ExperimentalAverageBelowCutoff, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_ModelAverageAboveCutoff, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_ModelAverageBelowCutoff, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_PositivesAboveThreshold, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_NumberActives, OSTREAM, INDENT) << '\n';

      // return the stream
      return OSTREAM;
    }

    //! @brief evaluates whether the approximation can continue
    //! @return true, if the approximation can continue - otherwise false
    bool ApproximatorAlignCutoff::CanContinue() const
    {
      return m_InternalIterate->CanContinue();
    }

    //! @brief determine the adjusted cutoff that maximizes the precision on the classification model
    void ApproximatorAlignCutoff::DetermineAdjustedCutoff()
    {
      const FeatureDataSet< float> predicted
      (
        m_LastInternallyGeneratedModel->PredictWithoutRescaling( *( m_InternalIterate->GetTrainingData()->GetFeaturesPtr()))
      );
      const FeatureDataSetInterface< float> &experimental( *( m_InternalIterate->GetTrainingData()->GetResultsPtr()));

      linal::MatrixConstReference< float> predicted_matrix( predicted.GetMatrix());

      BCL_Assert
      (
        experimental.GetNumberFeatures() == predicted.GetNumberFeatures(),
        "Experimental and predicted values do not have the same number of elements!"
      );

      BCL_Assert
      (
        experimental.GetFeatureSize() == predicted.GetFeatureSize(),
        "Experimental and predicted values do not have the same number of elements!"
      );

      // number of experimental values
      const size_t data_set_size( experimental.GetNumberFeatures());

      // number of predicted result columns
      const size_t result_size( predicted.GetFeatureSize());

      if( data_set_size == 0 || result_size == 0)
      {
        // no data available, no result column given
        m_ModelAverageBelowCutoff =  m_ExperimentalAverageBelowCutoff - m_ExperimentalCutoff + m_ModelCutoff;
        m_ModelAverageAboveCutoff = m_ExperimentalAverageAboveCutoff - m_ExperimentalCutoff + m_ModelCutoff;
        return;
      }

      for( size_t result_number( 0); result_number < result_size; ++result_number)
      {
        float best_precision( 0);

        // list of pairs with pred and exp values for ROC Curve
        storage::List< storage::Pair< double, double> > values_predicted_experimental;

        // iterate over all experimental and predicted values and check for accuracy
        for( size_t counter( 0); counter < data_set_size; ++counter)
        {
          values_predicted_experimental.PushBack
          (
            storage::Pair< double, double>
            (
              predicted( counter)( result_number),
              experimental( counter)( result_number)
            )
          );
        }

        // create roc curve according to cutoff
        math::ROCCurve roc_curve( values_predicted_experimental, m_ExperimentalCutoff, m_PositivesAboveThreshold);

        const storage::Vector< math::ROCCurve::Point> &sorted_counts( roc_curve.GetSortedCounts());
        size_t best_cutoff_count( 0);

        for
        (
          storage::Vector< math::ROCCurve::Point>::const_iterator
            itr( sorted_counts.Begin()), itr_end( sorted_counts.End());
          itr != itr_end;
          ++itr
        )
        {
          // # predicted (true & false) positives
          const size_t predicted( itr->GetNumberPredictedPositives());

          // require at least the actives are predicted before computing information gain ratio to avoid artifacts
          if( predicted < m_NumberActives)
          {
            continue;
          }

          math::ContingencyMatrix matrix( itr->GetContingencyMatrix( sorted_counts.LastElement()));

          // determine best seen precision value so far
          const float precision( matrix.GetInformationGainRatio());
          if( precision > best_precision)
          {
            best_precision = precision;
            best_cutoff_count = predicted;
            // set best cutoff seen so far as model cutoff
            m_ModelCutoff = itr->GetCutoff();
          }
        }

        // compute the average values above and below the cutoff for the model's values
        math::RunningAverage< float> average_above_cutoff, average_below_cutoff;
        for( size_t result_row_number( 0); result_row_number < data_set_size; ++result_row_number)
        {
          const float value( predicted_matrix( result_row_number, result_number));
          if( value > m_ModelCutoff)
          {
            average_above_cutoff += value;
          }
          else
          {
            average_below_cutoff += value;
          }
        }

        m_ModelAverageAboveCutoff = average_above_cutoff;
        m_ModelAverageBelowCutoff = average_below_cutoff;

        //
        if( average_above_cutoff.GetWeight() < 0.5)
        {
          // nothing above the cutoff was found
          m_ModelAverageAboveCutoff = m_ExperimentalAverageAboveCutoff - m_ExperimentalCutoff + m_ModelCutoff;
        }
        else if( average_below_cutoff.GetWeight() < 0.5)
        {
          m_ModelAverageBelowCutoff = m_ExperimentalAverageBelowCutoff - m_ExperimentalCutoff + m_ModelCutoff;
        }

        BCL_MessageStd
        (
          "AdjustedCutoff: " + util::Format()( m_ModelCutoff)
          + " best precision:" + util::Format()( best_precision)
          + " at fraction: " + util::Format()( float( best_cutoff_count) / float( data_set_size))
          + " av below: " + util::Format()( m_ModelAverageBelowCutoff)
          + " av above: " + util::Format()( m_ModelAverageAboveCutoff)
          + " # above: " + util::Format()( average_above_cutoff.GetWeight())
          + " # below: " + util::Format()( average_below_cutoff.GetWeight())
        );
      }

      // if the model does not predict any actives, then set the model slope above ave equal to the experimental slope
      if( math::EqualWithinTolerance( m_ModelAverageAboveCutoff, m_ModelCutoff, GetNumericalTolerance()))
      {
        m_ModelAverageAboveCutoff = m_ModelCutoff + m_ExperimentalAverageAboveCutoff - m_ExperimentalCutoff;
      }

      if( math::EqualWithinTolerance( m_ModelAverageBelowCutoff, m_ModelCutoff, GetNumericalTolerance()))
      {
        m_ModelAverageBelowCutoff = m_ModelCutoff + m_ExperimentalAverageBelowCutoff - m_ExperimentalCutoff;
      }
    }

  } // namespace model
} // namespace bcl
