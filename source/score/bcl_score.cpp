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
#include "score/bcl_score.h"

// includes from bcl - sorted alphabetically
#include "bcl_version.h"
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_command_state.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter.h"
#include "command/bcl_command_parameter_check_file_in_search_path.h"
#include "util/bcl_util_runtime_environment_interface.h"

// external includes - sorted alphabetically

//path for histograms to calculate energy functions
#if defined (__MINGW32__)
  #define BCL_HISTOGRAM_PATH "histogram/"
#elif defined (__GNUC__)
  #define BCL_HISTOGRAM_PATH "/dors/meilerlab/apps/bcl/histograms/rev_5265/"
#elif defined (_MSC_VER)
  #define BCL_HISTOGRAM_PATH "../../../histogram/"
#endif

namespace bcl
{
  namespace score
  {

    //! flag to change histogram path
    util::ShPtr< command::FlagInterface> &GetHistogramPathFlag()
    {
      static util::ShPtr< command::FlagInterface> s_histogram_path_flag
      (
        new command::FlagStatic
        (
          "histogram_path",
          "change path for reading and writing histograms for scores",
          command::Parameter
          (
            "path",
            "relative or absolute path",
            command::ParameterCheckFileInSearchPath
            (
              "histogram",
              Score::GetDefaultHistogramPath(),
              io::Directory::e_Dir
            ),
            ""
          )
        )
      );

      return s_histogram_path_flag;
    }

    //! @brief get the default histogram path
    //! @return the default histogram path - depending on release vs license vs normal
    const std::string &Score::GetDefaultHistogramPath()
    {
      static const std::string s_default_histogram_path
      (
        GetVersion().IsLicense() ? GetVersion().GetInstallPrefix() + "/histogram/" : std::string( BCL_HISTOGRAM_PATH)
      );

      return s_default_histogram_path;
    }

    //! @brief enum, that adds histogram params flag to default app flags
    static const util::ShPtr< command::FlagInterface> e_HistogramPathFlag
    (
      command::GetAppDefaultFlags().AddDefaultFlag( GetHistogramPathFlag(), command::e_Score)
    );

    //! @brief identifier for the name space
    //! @return the name of the namespace
    const std::string &GetNamespaceIdentifier()
    {
      static const std::string *s_namespace_name( new std::string( util::ExtractNamespaceIdentifier( __PRETTY_FUNCTION__)));
      return *s_namespace_name;
    }

    //! @brief given a FILENAME, the histogram path is prepended to the filename
    //! @param FILENAME the filename to a histogram that is used in one of the scores
    //! @return string with histogrampath/FILENAME
    std::string Score::AddHistogramPath( const std::string &FILENAME)
    {
      if( command::CommandState::IsInStaticInitialization())
      {
        return "histogram/" + FILENAME;
      }
      static const std::string s_hist_dir
      (
        GetHistogramPathFlag()->GetFlag()
        ? GetHistogramPathFlag()->GetFirstParameter()->GetValue()
        : command::ParameterCheckFileInSearchPath
        (
          "histogram",
          Score::GetDefaultHistogramPath(),
          io::Directory::e_Dir
        ).FindFile( "")
      );
      const std::string resolved_filename
      (
        util::GetRuntimeEnvironment().ResolveFileName( s_hist_dir + FILENAME)
      );

      BCL_Assert
      (
        !resolved_filename.empty(),
        "unable to resolve filename: " + GetHistogramPathFlag()->GetFirstParameter()->GetValue() + FILENAME
      );

      return resolved_filename;
    }

  } // namespace score
} // namespace bcl
