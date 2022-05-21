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
#include "model/bcl_model_approximator_base.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    //! @detail sets the improvement type of the tracker to SmallerIsBetter
    ApproximatorBase::ApproximatorBase() :
      opti::ApproximatorModularBase< util::ShPtr< Interface>, float>( opti::s_NumberImprovementTypes),
      m_TrainingContinued( false),
      m_ObjectiveFunction( new ObjectiveFunctionWrapper)
    {
    }

    //! @brief construct from objective function
    //! @detail sets the improvement type of the tracker to SmallerIsBetter
    //! @param OBJECTIVE_FUNCTION objective function used to evaluate a monitoring data set
    ApproximatorBase::ApproximatorBase( const util::ShPtr< ObjectiveFunctionWrapper> &OBJECTIVE_FUNCTION) :
      opti::ApproximatorModularBase< util::ShPtr< Interface>, float>
      (
        OBJECTIVE_FUNCTION.IsDefined()
        ? OBJECTIVE_FUNCTION->GetImprovementType()
        : opti::s_NumberImprovementTypes
      ),
      m_TrainingContinued( false),
      m_ObjectiveFunction
      (
        OBJECTIVE_FUNCTION.IsDefined()
        ? OBJECTIVE_FUNCTION
        : util::ShPtr< ObjectiveFunctionWrapper>( new ObjectiveFunctionWrapper)
      )
    {
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns a shared pointer to the training data set
    //! @return shared pointer to the data set interface with training data
    const util::ShPtr< descriptor::Dataset> &ApproximatorBase::GetTrainingData() const
    {
      return m_TrainingData;
    }

    //! @brief get objective function to evaluate a monitoring data set
    //! @return data set interface with training data
    const util::ShPtr< ObjectiveFunctionWrapper> &ApproximatorBase::GetObjectiveFunction() const
    {
      return m_ObjectiveFunction;
    }

    //! @brief sets the objective function used to evaluate a monitoring data set
    //! @param OBJECTIVE_FUNCTION objective function used to evaluate a monitoring data set
    void ApproximatorBase::SetObjectiveFunction( const util::ShPtr< ObjectiveFunctionWrapper> &OBJECTIVE_FUNCTION)
    {
      m_ObjectiveFunction = OBJECTIVE_FUNCTION;
      this->GetTracker() = opti::Tracker< util::ShPtr< Interface>, float>( m_ObjectiveFunction->GetImprovementType());
    }

    //! @brief returns the flag for continuation of training
    //! @return flag for continuation of training
    bool ApproximatorBase::IsTrainingContinued() const
    {
      return m_TrainingContinued;
    }

    //! @brief sets the flag for continuation of training
    //! @param flag for continuation of training
    void ApproximatorBase::SetTrainingContinued( const bool &CONTINUED_TRAINING)
    {
      m_TrainingContinued = CONTINUED_TRAINING;
    }

    //! @brief returns a shared pointer to the rescale function to transfer data sets into a scaled feature space
    //! @return shared pointer to rescale function to rescale into normalized feature space
    util::ShPtr< RescaleFeatureDataSet> ApproximatorBase::GetRescaleFeatureDataSet() const
    {
      return
        m_TrainingData.IsDefined()
        ? m_TrainingData->GetFeaturesPtr()->GetScaling()
        : util::ShPtr< RescaleFeatureDataSet>();
    }

    //! @brief returns a shared pointer to the rescale function to transfer data sets into a scaled result space
    //! @return shared pointer to the rescale function to rescale into normalized result space
    util::ShPtr< RescaleFeatureDataSet> ApproximatorBase::GetRescaleResultDataSet() const
    {
      return
        m_TrainingData.IsDefined()
        ? m_TrainingData->GetResultsPtr()->GetScaling()
        : util::ShPtr< RescaleFeatureDataSet>();
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERR_STREAM stream to write errors out to
    bool ApproximatorBase::ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERR_STREAM)
    {
      if( m_ObjectiveFunction.IsDefined())
      {
        this->GetTracker() = opti::Tracker< util::ShPtr< Interface>, float>( m_ObjectiveFunction->GetImprovementType());
      }
      return true;
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ApproximatorBase::Read( std::istream &ISTREAM)
    {
      // call Read of base class
      opti::ApproximatorModularBase< util::ShPtr< Interface>, float>::Read( ISTREAM);

      // All members of ApproximatorBase are intended for being reset whenever the object is read back in for
      // continuation, since the purpose of reading it back in may be to change the training data or objective function
      m_TrainingContinued = true;

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &ApproximatorBase::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // call Write of base class
      opti::ApproximatorModularBase< util::ShPtr< Interface>, float>::Write( OSTREAM, INDENT) << '\n';

      // All members of ApproximatorBase are intended for being reset whenever the object is read back in for
      // continuation, since the purpose of reading it back in may be to change the training data or objective function

      // end
      return OSTREAM;
    }

  } // namespace model
} // namespace bcl
