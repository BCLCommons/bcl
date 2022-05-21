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

// includes from bcl - sorted alphabetically
#include "app/bcl_app_apps.h"
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_command.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter_check_serializable.h"
#include "io/bcl_io_file.h"
#include "model/bcl_model_feature_similarity_interface.h"
#include "model/bcl_model_retrieve_data_set_base.h"
#include "util/bcl_util_implementation.h"
#include "util/bcl_util_stopwatch.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace app
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class DescriptorDatasetSimilarityMeasures
    //! @brief generates similarity measures of small molecules based on feature vectors or of feature vectors themselves
    //! @details similarity measures currently implemented are Tanimoto, Manhattan and Euclidean distance, Cosine,
    //!          and Dice
    //!
    //! @author loweew
    //! @date Sep 10, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API DescriptorDatasetSimilarityMeasures :
      public InterfaceRelease
    {

    private:

    //////////
    // data //
    //////////

      //! flag for input dataset
      util::ShPtr< command::FlagInterface> m_InputDataSet;

      //! flag for output table file name
      util::ShPtr< command::FlagInterface> m_OutputFilenameFlag;

      //! flag for code object file for features vector generation
      util::ShPtr< command::FlagInterface> m_InputCodeFlag;
      util::ShPtr< command::FlagInterface> m_OutputCodeFlag;

      //! flag for the similarity measure label i.e. 'DatasetSimilarityMeasures(measurement=Tanimoto)'
      util::ShPtr< command::FlagInterface> m_MeasureFlag;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      DescriptorDatasetSimilarityMeasures();

      //! @brief Clone function
      //! @return pointer to new DescriptorDatasetSimilarityMeasures
      DescriptorDatasetSimilarityMeasures *Clone() const
      {
        return new DescriptorDatasetSimilarityMeasures( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief get a description for the app
      //! @return a brief (no more than 3 line) description for the application
      std::string GetDescription() const
      {
        return "Generate similarity measures of small molecule feature sets.";
      }

      //! @brief returns readme information
      //! @return string containing information about application
      const std::string &GetReadMe() const
      {
        static std::string s_read_me =
            "DatasetSimilarityMeasures is an application that generates similarity measures of small molecules "
            "based on feature vectors or of feature vectors themselves. Similarity measures currently implemented are: "
            "Tanimoto, Manhattan and Euclidean distance, Cosine, and Dice";
        return s_read_me;
      }

      //! @brief initializes the command object for that executable
      util::ShPtr< command::Command> InitializeCommand() const
      {
        // initialize command to be returned
        util::ShPtr< command::Command> sp_cmd( new command::Command());

        // insert all the flags and params
        sp_cmd->AddFlag( m_InputDataSet);
        sp_cmd->AddFlag( m_OutputFilenameFlag);
        sp_cmd->AddFlag( m_InputCodeFlag);
        sp_cmd->AddFlag( m_OutputCodeFlag);
        sp_cmd->AddFlag( m_MeasureFlag);

        // default flags
        command::GetAppDefaultFlags().AddDefaultCommandlineFlags( *sp_cmd);

        // return assembled Command object
        return sp_cmd;
      }

      //! @brief the Main function
      //! @return error code - 0 for success
      int Main() const
      {
        // initialize implementations
        util::Implementation< model::RetrieveDataSetBase> retriever( m_InputDataSet->GetFirstParameter()->GetValue());

        // read in feature / result codes
        retriever->SelectFeaturesGivenFilenameFlag( *m_InputCodeFlag);
        retriever->SelectResultsGivenFilenameFlag( *m_OutputCodeFlag);

        // create the input matrix
        util::ShPtr< descriptor::Dataset> dataset( retriever->GenerateDataSet());
        linal::MatrixConstReference< float> input_matrix( dataset->GetFeaturesReference());
        // check whether dataset is not empty
        BCL_Assert( input_matrix.GetNumberOfElements(), "dataset was empty!");

        util::Implementation< model::FeatureSimilarityMeasuresInterface< float> > measures
        (
          m_MeasureFlag->GetFirstParameter()->GetValue()
        );

        // allocate host memory for result
        linal::Matrix< float> results( input_matrix.GetNumberRows(), input_matrix.GetNumberRows());
        util::Stopwatch timer;
        // if the data sets won't fit on the GPU, split
        timer.Reset();
        timer.Start();

        // perform operation
        results = measures->operator()( input_matrix);

        // write out result
        io::OFStream out;
        io::File::MustOpenOFStream( out, m_OutputFilenameFlag->GetFirstParameter()->GetValue());
        out << results;
        io::File::CloseClearFStream( out);

        // end
        return 0;
      }

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

      // instantiate enumerator for DescriptorDatasetSimilarityMeasures class
      static const ApplicationType DescriptorDatasetSimilarityMeasures_Instance;

    }; // class DescriptorDatasetSimilarityMeasures

    //! @brief standard constructor
    DescriptorDatasetSimilarityMeasures::DescriptorDatasetSimilarityMeasures() :
      m_InputDataSet
      (
        new command::FlagStatic
        (
          "input",
          "flag for input data set",
          command::Parameter
          (
            "input",
            "",
            command::ParameterCheckSerializable( util::Implementation< model::RetrieveDataSetBase>())
          )
        )
      ),
      m_OutputFilenameFlag
      (
        new command::FlagStatic
        (
          "output_filename", "flag for output file name",
          command::Parameter
          (
            "output_filename", "flag for output file name", "output.mat"
          )
        )
      ),
      m_InputCodeFlag
      (
        new command::FlagStatic
        (
          "input_code", "flag for code object file name",
          command::Parameter
          (
            "input_code", "flag for code object file name", ""
          )
        )
      ),
      m_OutputCodeFlag
      (
        new command::FlagStatic
        (
          "output_code", "flag for code object file name",
          command::Parameter
          (
            "output_code", "flag for code object file name", ""
          )
        )
      ),
      m_MeasureFlag
      (
        new command::FlagStatic
        (
          "measure", "flag for the measure to be calculated (Tanimoto, Dice, Cosine, Euclidean, Manhattan)",
          command::Parameter
          (
            "measure",
            "flag for the measure to be calculated (Tanimoto, Dice, Cosine, Euclidean, Manhattan)",
            command::ParameterCheckSerializable( util::Implementation< model::FeatureSimilarityMeasuresInterface< float> >()),
            ""
          )
        )
      )
    {
    }

    const ApplicationType DescriptorDatasetSimilarityMeasures::DescriptorDatasetSimilarityMeasures_Instance
    (
      GetAppGroups().AddAppToGroup( new DescriptorDatasetSimilarityMeasures(), GetAppGroups().e_Descriptor)
    );

  } // namespace app
} // namespace bcl
