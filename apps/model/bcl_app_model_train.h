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

#ifndef BCL_APP_MODEL_TRAIN_H_
#define BCL_APP_MODEL_TRAIN_H_

// include the namespace header
#include "app/bcl_app_apps.h"

// include forward headers
#include "descriptor/bcl_descriptor.fwd.hh"
#include "model/bcl_model.fwd.hh"
#include "opti/bcl_opti.fwd.hh"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_flag_interface.h"

namespace bcl
{
  namespace app
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ModelTrain
    //!
    //! @brief Application for training a model (model::Interface e.g., NeuralNetwork, SVM). Using the approximator
    //!        framework this application reads in necessary object files, as well as a container of dataset balanced
    //!        representing all the chunks involed in a cross-validation training run. Given the ids for independent
    //!        and monitoring data chunk one iteration in a cross-validation scheme will be executed.
    //!
    //! @see @link example_app_model_train.cpp @endlink
    //!
    //! @author butkiem1, mendenjl
    //! @date 12/17/2009
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ModelTrain :
      public InterfaceRelease
    {
    private:

    //////////
    // data //
    //////////

      //! iterate label
      util::ShPtr< command::ParameterInterface> m_Iterator;

      //! whether to suppress output objective function results while training
      util::ShPtr< command::FlagInterface> m_SuppressProgressOutput;
      //! maximum time, in minutes; the training is stopped at the next iteration thereafter
      util::ShPtr< command::FlagInterface> m_FlagTerminateAfterTimeMinutes;
      //! maximum # of iterations; a termination criteria that is added to the terminate file
      util::ShPtr< command::FlagInterface> m_FlagMaxNumberIterations;
      //! maximum # of iterations since the last improved step
      util::ShPtr< command::FlagInterface> m_FlagMaxIterationsWithoutImprovement;
      //! Number of rounds over which to take the (triangularly-weighted) average result to compute the current result
      //! This is useful for noisy objective functions (e.g. Enrichment), to avoid spurious improvements being accepted
      util::ShPtr< command::FlagInterface> m_FlagTrackerResultAverageRounds;

      //! objective function label for final evaluation
      util::ShPtr< command::FlagInterface> m_FlagFinalObjectiveFunction;

      //! flag for writing out an iterate instance for continued training
      util::ShPtr< command::FlagInterface> m_FlagContinuedTraining;
      //! flag for reading in an iterate instance for continued training
      util::ShPtr< command::FlagInterface> m_FlagReadApproximator;

      //! flag for the training data set
      util::ShPtr< command::FlagInterface> m_FlagTrainingDataSet;
      //! flag for the monitoring data set
      util::ShPtr< command::FlagInterface> m_FlagMonitoringDataSet;
      //! flag for the independent data set
      util::ShPtr< command::FlagInterface> m_FlagIndependentDataSet;
      //! flag for the list of feature names/types
      util::ShPtr< command::FlagInterface> m_FlagFeatureCode;
      //! flag for the list of result names/types
      util::ShPtr< command::FlagInterface> m_FlagResultCode;
      //! flag for the list of id names/types
      util::ShPtr< command::FlagInterface> m_FlagIdCode;
      //! flag for whether to print the predicted values for training
      util::ShPtr< command::FlagInterface> m_FlagPrintTrainingPredictions;
      //! flag for whether to print the predicted values for training
      util::ShPtr< command::FlagInterface> m_FlagPrintMonitoringPredictions;
      //! flag for whether to print the predicted values for training
      util::ShPtr< command::FlagInterface> m_FlagPrintIndependentPredictions;
      //! flag for storing the final trained model
      util::ShPtr< command::FlagInterface> m_FlagModelInterfaceStorage;
      //! flag for storing meta data information for descriptor selection
      util::ShPtr< command::FlagInterface> m_FlagMetaDataStorage;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! default constructor
      ModelTrain();

    public:

      //! @brief Clone function
      //! @return pointer to new ModelTrain instance
      ModelTrain *Clone() const
      {
        return new ModelTrain( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      // This function is used to retrieve prior names for applications that are still valid (for backwards compatibility)
      // but which should not be displayed, e.g. for help
      storage::Vector< std::string> GetDeprecatedAppNames() const;

      //! @brief get a description for the app
      //! @return a brief (no more than 3 line) description for the application
      std::string GetDescription() const;

      //! @brief returns readme information
      //! @return string containing information about application
      const std::string &GetReadMe() const;

      //! @brief returns web text information
      //! @return text (html allowed but not required) that will be displayed on the website
      //! @note use confluence syntax, e.g. !my_image.png! to refer to an image under documentation/image
      const std::string &GetWebText() const;

      //! @brief initializes the command object for that executable
      util::ShPtr< command::Command> InitializeCommand() const;

    public:

      //! @brief the Main function
      //! @return error code - 0 for success
      int Main() const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM)
      {
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        return OSTREAM;
      }

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      //! @brief Create termination criterion from command line flag
      opti::CriterionCombine< util::ShPtr< model::Interface>, float> CreateTerminationCriterion() const;

    public:

      // instantiate enumerator for ModelTrain class
      static const ApplicationType ModelTrain_Instance;

    }; // ModelTrain

  } // namespace app
} // namespace bcl

#endif // BCL_APP_MODEL_TRAIN_H_
