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

#ifndef BCL_APP_MODEL_PREDICTION_MERGE_H_
#define BCL_APP_MODEL_PREDICTION_MERGE_H_

// include the namespace header

// include forward headers
#include "linal/bcl_linal.fwd.hh"
#include "model/bcl_model.fwd.hh"

// includes from bcl - sorted alphabetically
#include "app/bcl_app_apps.h"
#include "command/bcl_command_flag_interface.h"

namespace bcl
{
  namespace app
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ModelPredictionMerge
    //! @brief ModelPredictionMerge handles prediction outputs written out by app::TrainModel
    //!
    //! @details This application merges or appends prediction matrices to create a cross-validation matrix
    //!          eg. 2 column matrix with experimental and predicted values
    //!
    //! @author butkiem1
    //!
    //! @date 06/11/2010
    //!
    //! @see @link example_app_model_prediction_merge.cpp @endlink
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ModelPredictionMerge :
      public InterfaceRelease
    {

    private:

    //////////
    // data //
    //////////

      //! filenames with input prediction matrices
      util::ShPtr< command::FlagInterface> m_InputFilenames;

      //! flag to append predictions (usually from different datasets) rather than averaging them
      util::ShPtr< command::FlagInterface> m_AppendMerge;

      //! flag to use a model storage to compute all the input filenames, which to merge, which to append, etc., rather
      //! than m_InputFilenames and m_AppendMerge
      util::ShPtr< command::FlagInterface> m_ModelStorage;

      //! filename for output matrix
      util::ShPtr< command::FlagInterface> m_OutputFilename;

      //! flag to take median prediction (rather than mean)
      util::ShPtr< command::FlagInterface> m_Median;

      //! flag to take jury prediction (rather than mean)
      mutable util::ShPtr< command::FlagInterface> m_Jury;

      //! flag to take min prediction (rather than mean)
      util::ShPtr< command::FlagInterface> m_Min;

      //! flag to take max prediction (rather than mean)
      util::ShPtr< command::FlagInterface> m_Max;

      //! flag to take mean of LocalPPV values.
      util::ShPtr< command::FlagInterface> m_LocalPPV;

    ///////////////////////////////////
    // construction and destruction //
    ///////////////////////////////////

      //! default constructor
      ModelPredictionMerge();

    public:

      //! @brief Clone function
      //! @return pointer to new FoldProtein
      ModelPredictionMerge *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

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
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief compute median for each result column, output into a matrix
      //! @param RESULTS vector of predictions, all corresponding to the same dataset
      //! @return median prediction values for each result column
      linal::Matrix< float> ComputeMedian( const storage::Vector< linal::Matrix< float> > &RESULTS) const;

      //! @brief compute jury for each result column, output into a matrix
      //! @param PREDICTIONS vector of predictions, all corresponding to the same dataset
      //! @param RESULTS actual results
      //! @param RESULT_IDS ids for each result, needed by some objective functions
      //! @return jury prediction values for each result column
      linal::Matrix< float> ComputeJury
      (
        const storage::Vector< linal::Matrix< float> > &PREDICTIONS,
        const model::FeatureDataSet< float> &RESULTS,
        const model::FeatureDataSet< char> &RESULT_IDS
      ) const;

    public:

      // instantiate enumerator for PrepareSmallMoleculeEnsemble class
      static const ApplicationType ModelPredictionMerge_Instance;

    }; // ModelPredictionMerge

  } // namespace app
} // namespace bcl

#endif // BCL_APP_MODEL_PREDICTION_MERGE_H_
