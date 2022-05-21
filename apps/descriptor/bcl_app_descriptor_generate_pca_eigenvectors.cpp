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
#include "command/bcl_command_flag_dynamic.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter_check_file_existence.h"
#include "command/bcl_command_parameter_check_serializable.h"
#include "io/bcl_io_file.h"
#include "linal/bcl_linal_principal_component_analysis.h"
#include "math/bcl_math_running_average.h"
#include "math/bcl_math_template_instantiations.h"
#include "model/bcl_model_objective_function_wrapper.h"
#include "model/bcl_model_retrieve_data_set_base.h"
#include "storage/bcl_storage_template_instantiations.h"
#include "util/bcl_util_stopwatch.h"

namespace bcl
{
  namespace app
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class DescriptorGeneratePCAEigenVectors
    //!
    //! @brief this application generates the eigenvectors matrix need for extracting principle components from data set
    //!
    //! @author loweew
    //! @date 06/29/2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API DescriptorGeneratePCAEigenVectors :
      public InterfaceRelease
    {
    private:

    //////////
    // data //
    //////////

      //! flag for the training data set
      util::ShPtr< command::FlagInterface> m_FlagTrainingDataSet;
      //! flag for the list of feature names/types
      util::ShPtr< command::FlagInterface> m_FlagFeatureCode;
      //! flag for the list of result names/types
      util::ShPtr< command::FlagInterface> m_FlagOutputFilename;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! default constructor
      DescriptorGeneratePCAEigenVectors();

    public:

      //! @brief Clone function
      //! @return pointer to new FoldProtein
      DescriptorGeneratePCAEigenVectors *Clone() const
      {
        return new DescriptorGeneratePCAEigenVectors( *this);
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

      //! @brief get a description for the app
      //! @return a brief (no more than 3 line) description for the application
      std::string GetDescription() const
      {
        return "Generate descriptor set eigenvectors";
      }

      //! @brief returns readme information
      //! @return string containing information about application
      const std::string &GetReadMe() const
      {
        static std::string s_read_me =
            "GeneratePCAEigenVectors is an application that generates the "
            "eigenvectors matrix needed for extracting principle components from "
            "descriptor datasets.";
        return s_read_me;
      }

      //! @brief initializes the command object for that executable
      util::ShPtr< command::Command> InitializeCommand() const
      {
        util::ShPtr< command::Command> sp_cmd( new command::Command());

        // parameters

        // training, monitoring, independent
        sp_cmd->AddFlag( m_FlagTrainingDataSet);

        // code object files
        sp_cmd->AddFlag( m_FlagFeatureCode);
        sp_cmd->AddFlag( m_FlagOutputFilename);

        // add default bcl parameters
        command::GetAppDefaultFlags().AddDefaultCommandlineFlags( *sp_cmd);

        // return assembled Command object
        return sp_cmd;
      } // InitializeCommand

    ////////////////////
    // helper methods //
    ////////////////////

    private:

    public:

      //! @brief the Main function
      //! @return error code - 0 for success
      int Main() const
      {
        // stopwatch
        util::Stopwatch timer;
        timer.Start();

        // initialize monitor data
        util::Implementation< model::RetrieveDataSetBase> dataset_retriever
        (
          util::ObjectDataLabel( m_FlagTrainingDataSet->GetFirstParameter()->GetValue())
        );
        // set up the feature result codes
        dataset_retriever->SelectFeaturesGivenFilenameFlag( *m_FlagFeatureCode);

        // create sh-ptrs to the different datasets
        linal::Matrix< float> features( dataset_retriever->GenerateDataSet()->GetFeaturesReference());

        // check whether monitoring dataset is not empty
        BCL_Assert( features.GetNumberOfElements(), "dataset is empty!");

        model::RescaleFeatureDataSet rescale_in
        (
          features,
          math::Range< float>( -1.0f, 1.0f),
          model::RescaleFeatureDataSet::e_AveStd
        );

        rescale_in.RescaleMatrix( features);

        BCL_MessageStd
        (
          "Dataset retrieved and rescaled in " + util::Format()( timer.GetTotalTime().GetSeconds())
          + " [sec]\nNumber of features: " + util::Format()( features.GetNumberRows())
          + "\nNumber of descriptor values: " + util::Format()( features.GetNumberCols())
        );

        timer.Reset();
        timer.Start();

        storage::Pair< linal::Matrix< float>, linal::Vector< float> > eig_vecs_and_values
        (
          linal::PrincipalComponentAnalysis< float>::GetSortedEigenVectorsValues( features)
        );

        // READ IN
        BCL_MessageStd
        (
          "total time running pca: " + util::Format()( timer.GetTotalTime().GetSeconds()) + " [sec]"
        );

        io::OFStream output;
        io::File::MustOpenOFStream( output, m_FlagOutputFilename->GetFirstParameter()->GetValue());
        io::Serialize::Write( rescale_in, output);
        io::Serialize::Write( eig_vecs_and_values, output);
        io::File::CloseClearFStream( output);

        linal::Matrix< float> abs_eigenvec( eig_vecs_and_values.First());
        const linal::Vector< float> &eigvenval( eig_vecs_and_values.Second());
        math::Absolute( abs_eigenvec);
        math::RunningAverage< linal::Vector< float> > ave_col_importance;
        for( size_t i( 0), n_rows( abs_eigenvec.GetNumberRows()); i < n_rows; ++i)
        {
          ave_col_importance.AddWeightedObservation( abs_eigenvec.GetRow( i), eigvenval( i));
        }

        model::FeatureLabelSet labels_split( dataset_retriever->GetFeatureLabelsWithSizes().SplitFeatureLabelSet( true));
        io::File::MustOpenOFStream( output, m_FlagOutputFilename->GetFirstParameter()->GetValue() + ".score");
        for( size_t i( 0), n_rows( abs_eigenvec.GetNumberRows()); i < n_rows; ++i)
        {
          output << ave_col_importance.GetAverage()( i) << '\t' << labels_split.GetMemberLabels()( i).ToString() << '\n';
        }
        io::File::CloseClearFStream( output);
        // end
        return 0;
      } // Main

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

      // instantiate enumerator for PrepareSmallMoleculeEnsemble class
      static const ApplicationType DescriptorGeneratePCAEigenVectors_Instance;

    }; // GeneratePCAEigenVectors

    //! @brief standard constructor
    DescriptorGeneratePCAEigenVectors::DescriptorGeneratePCAEigenVectors() :
      m_FlagTrainingDataSet
      (
        new command::FlagStatic
        (
          "training",
          "method to retrieve the training dataset",
          command::Parameter
          (
            "",
            "",
            command::ParameterCheckSerializable( util::Implementation< model::RetrieveDataSetBase>())
          )
        )
      ),
      m_FlagFeatureCode
      (
        new command::FlagDynamic
        (
          "feature_labels",
          "file containing the feature label",
          command::Parameter
          (
            "",
            "",
            command::ParameterCheckFileExistence()
          ),
          0,
          1
        )
      ),
      m_FlagOutputFilename
      (
        new command::FlagStatic
        (
          "output_filename", "the output file name",
          command::Parameter
          (
            "output_filename", "output file name"
          )
        )
      )
    {
    }

    const ApplicationType DescriptorGeneratePCAEigenVectors::DescriptorGeneratePCAEigenVectors_Instance
    (
      GetAppGroups().AddAppToGroup( new DescriptorGeneratePCAEigenVectors(), GetAppGroups().e_Descriptor)
    );

  } // namespace app
} // namespace bcl
