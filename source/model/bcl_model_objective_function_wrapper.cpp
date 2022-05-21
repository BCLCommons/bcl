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
#include "model/bcl_model_objective_function_wrapper.h"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_running_average_sd.h"
#include "model/bcl_model_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> ObjectiveFunctionWrapper::s_Instance
    (
      GetObjectInstances().AddInstance( new ObjectiveFunctionWrapper())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    ObjectiveFunctionWrapper::ObjectiveFunctionWrapper() :
      m_Data(),
      m_Specialization()
    {
    }

    //! @brief constructor with parameters
    //! @param SPECIALIZATION specify objective function implementation
    ObjectiveFunctionWrapper::ObjectiveFunctionWrapper
    (
      const util::Implementation< ObjectiveFunctionInterface> &SPECIALIZATION
    ) :
      m_Data(),
      m_Specialization( SPECIALIZATION)
    {
    }

    //! @brief constructor with parameters
    //! @param DATA data set that is used for testing the progress of the model::Interface
    //! @param SPECIALIZATION specify objective function implementation
    ObjectiveFunctionWrapper::ObjectiveFunctionWrapper
    (
      util::ShPtr< descriptor::Dataset> &DATA,
      const util::Implementation< ObjectiveFunctionInterface> &SPECIALIZATION
    ) :
      m_Data( DATA),
      m_Specialization( SPECIALIZATION)
    {
    }

    //! @brief Clone function
    //! @return pointer to new ObjectiveFunctionWrapper
    ObjectiveFunctionWrapper *ObjectiveFunctionWrapper::Clone() const
    {
      return new ObjectiveFunctionWrapper( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ObjectiveFunctionWrapper::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief set data and rescaling for support vector model
    //! @param DATA dataset of interest
    //! @param RESCALE_INPUT, RESCALE_OUTPUT the rescale functions
    void ObjectiveFunctionWrapper::SetData
    (
      util::ShPtr< descriptor::Dataset> &DATA,
      const util::ShPtr< RescaleFeatureDataSet> &RESCALE_INPUT
    )
    {
      BCL_Assert( DATA.IsDefined(), "Data set undefined!");
      m_Data = DATA;

      if( RESCALE_INPUT.IsDefined())
      {
        m_Data->GetFeatures().Rescale( *RESCALE_INPUT);
      }
      if( m_Specialization.IsDefined())
      {
        if( DATA->GetIdsPtr().IsDefined())
        {
          m_Specialization->SetData( *DATA->GetResultsPtr(), *DATA->GetIdsPtr());
        }
        else
        {
          m_Specialization->SetData( *DATA->GetResultsPtr());
        }
      }
    }

    //! @brief get monitoring data for model::Interface
    const util::ShPtr< descriptor::Dataset> &ObjectiveFunctionWrapper::GetData() const
    {
      return m_Data;
    }

    //! @brief get the overall goal of the objective function
    //! @return the goal of the objective function
    ObjectiveFunctionInterface::Goal ObjectiveFunctionWrapper::GetGoalType() const
    {
      // return the improvement type; regression by default if no objective function is defined
      return m_Specialization.IsDefined() ? m_Specialization->GetGoalType() : ObjectiveFunctionInterface::e_Regression;
    }

    //! @brief determine what sign of the derivative of this objective function indicates improvement
    //! @return the sign of the derivative of this objective function indicates improvement
    opti::ImprovementType ObjectiveFunctionWrapper::GetImprovementType() const
    {
      // return improvement types; undefined if none is available
      return
        m_Specialization.IsDefined()
        ? m_Specialization->GetImprovementType()
        : opti::s_NumberImprovementTypes;
    }

    //! @brief get implementation used by the objective function
    const util::Implementation< ObjectiveFunctionInterface> &
    ObjectiveFunctionWrapper::GetImplementation() const
    {
      return m_Specialization;
    }

    //! @brief get the threshold, for classification type objectives
    //! @return the threshold, for classification type objectives
    float ObjectiveFunctionWrapper::GetThreshold() const
    {
      return m_Specialization->GetThreshold();
    }

    //! @brief get the parity, for rank classification type objectives
    //! @return the parity, for rank classification type objectives
    //! @note this has meaning only for rank-classification objectives; true means the objective is most interested
    //!       in prediction of values higher than the threshold, false means below threshold
    bool ObjectiveFunctionWrapper::GetRankingParity() const
    {
      return m_Specialization->GetRankingParity();
    }

    //! @brief get the desired hit rate (e.g. fraction of positive predictions)
    //! @return the desired hit rate, for rank classification type objectives
    //! @note this has meaning only for rank-classification objectives; 0.01 means, for example, that only the top 1%
    //!       of values will be considered
    float ObjectiveFunctionWrapper::GetDesiredHitRate() const
    {
      return m_Specialization->GetDesiredHitRate();
    }

    //! @brief test whether one objective function value is better than another
    //! @param NEW_RESULT, OLD_RESULT two objective function results
    //! @return true if NEW_RESULT is better than OLD_RESULT
    bool ObjectiveFunctionWrapper::TestWhetherResultsImproved( const float &NEW_RESULT, const float &OLD_RESULT)
    {
      if( GetImprovementType() == opti::e_LargerEqualIsBetter)
      {
        return NEW_RESULT > OLD_RESULT;
      }

      return NEW_RESULT < OLD_RESULT;
    }

    //! @brief alter an output descaling function based on the objective function
    //! @param RESCALING current rescaling function
    //! @return the optimized rescaling function
    util::ShPtr< RescaleFeatureDataSet> ObjectiveFunctionWrapper::OptimizeRescalingFunction
    (
      const util::ShPtr< RescaleFeatureDataSet> &RESCALE,
      const FeatureDataSet< float> &PREDICTIONS
    ) const
    {
      // this function
      // computes a rescaling function as f(x) = (x-adjusted cutoff)*slope+EXPERIMENTAL_CUTOFF
      // where adjusted cutoff is the predicted value at which FRACTION have been predicted
      // determine slope

      // well defined hit rate, determine rescaling functions for every cutoff
      storage::Vector< math::Range< float> > old_ranges( RESCALE->GetRescaleRanges());

      const size_t dataset_size( PREDICTIONS.GetNumberFeatures());

      // compute the desired number of predicted positives

      // get a matrix reference
      linal::MatrixConstReference< float> predictions( PREDICTIONS.GetMatrix());
      util::ShPtrVector< math::FunctionInterfaceSerializable< float, float> > rescalers( PREDICTIONS.GetFeatureSize());
      for( size_t result( 0), result_size( PREDICTIONS.GetFeatureSize()); result < result_size; ++result)
      {
        math::RunningAverageSD< float> average_sd;
        for( size_t row( 0); row < dataset_size; ++row)
        {
          average_sd += predictions( row, result);
        }

        // get the original range width
        const float original_range_width( old_ranges( result).GetWidth());

        // get the predicted range width
        const float predicted_range_width( 2.0 * math::Absolute( average_sd.GetStandardDeviation()));

        // compute the slope
        float slope( 1.0);
        if( original_range_width && predicted_range_width)
        {
          slope = original_range_width / predicted_range_width;
        }

        // compute the offset
        const float offset( -slope * average_sd.GetAverage());

        // create the linear function to map values in the old range onto the new range
        rescalers( result) = util::CloneToShPtr( math::LinearFunction( slope, offset));
      }

      // create a new rescaler with the updated ranges
      return util::ShPtr< RescaleFeatureDataSet>( new RescaleFeatureDataSet( *RESCALE, rescalers));
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief compute error between predicted values from FEATURE and RESULT
    //! @param MODEL NeuralNetwork to be monitored
    //! @return error between predicted values from FEATURE and RESULT
    float ObjectiveFunctionWrapper::operator()( const util::PtrInterface< Interface> &MODEL) const
    {
      // predict a list of experimental and predicted values and evaluate objective function
      if( m_Specialization.IsDefined() && m_Data.IsDefined() && m_Data->GetResultsPtr().IsDefined())
      {
        return m_Specialization->operator()( *m_Data->GetResultsPtr(), Predict( MODEL));
      }
      // no objective function: just return 0
      return 0.0;
    }

    //! @brief compute feature dataset with predicted values
    //! @param MODEL shptr on model::Interface that will be monitored
    //! @return predicted values from FEATURE and RESULT based on given model
    FeatureDataSet< float> ObjectiveFunctionWrapper::Predict( const util::PtrInterface< Interface> &MODEL) const
    {
      // return list with predicted and experimental values
      return MODEL->PredictWithoutRescaling( m_Data->GetFeatures()).DeScale();
    }

    //! @brief compute feature dataset with predicted values
    //! @param PREDICTIONS predictions for which to evaluate this measure
    //!        PREDICTIONS may be rescaled, as necessary, for this operation, but it will be returned
    //!        with the same scaling it was given during input
    //! @return predicted values from FEATURE and RESULT based on given model
    float ObjectiveFunctionWrapper::Evaluate( FeatureDataSet< float> &PREDICTIONS) const
    {
      // return list with predicted and experimental values
      if( m_Specialization.IsDefined())
      {
        // handle rank-classification tasks, in which scaling does not matter
        if
        (
          GetGoalType() == ObjectiveFunctionInterface::e_RankClassification
          || PREDICTIONS.HasSameScaling( *m_Data->GetResultsPtr())
        )
        {
          return m_Specialization->operator()( *m_Data->GetResultsPtr(), PREDICTIONS);
        }
        // If predictions is scaled, it must be descaled first
        if( PREDICTIONS.IsRescaled())
        {
          util::ShPtr< RescaleFeatureDataSet> rescaling( PREDICTIONS.GetScaling());
          PREDICTIONS.DeScale();
          float result( Evaluate( PREDICTIONS));
          PREDICTIONS.Rescale( *rescaling);
          return result;
        }
        if( m_Data->GetResultsPtr()->IsRescaled())
        {
          util::ShPtr< RescaleFeatureDataSet> rescaling( m_Data->GetResultsPtr()->GetScaling());
          m_Data->GetResults().DeScale();
          float result( Evaluate( PREDICTIONS));
          m_Data->GetResults().Rescale( *rescaling);
          return result;
        }
      }
      return 0.0;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ObjectiveFunctionWrapper::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Specialization, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &ObjectiveFunctionWrapper::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Specialization, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

  } // namespace model

} // namespace bcl
