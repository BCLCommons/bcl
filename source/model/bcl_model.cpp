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
#include "model/bcl_model.h"

// includes from bcl - sorted alphabetically
#include "bcl_version.h"
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter.h"
#include "command/bcl_command_parameter_check_file_in_search_path.h"
#include "util/bcl_util_runtime_environment_interface.h"

// external includes - sorted alphabetically

//path for models
#if defined (__MINGW32__)
  #define BCL_MODEL_PATH "model/"
#elif defined (__GNUC__)
  #define BCL_MODEL_PATH "/dors/meilerlab/apps/bcl/model/rev_5051/"
#elif defined (_MSC_VER)
  #define BCL_MODEL_PATH "../../../model/"
#endif

namespace bcl
{
  namespace model
  {

    //! flag to change model path
    util::ShPtr< command::FlagInterface> &Model::GetModelPathFlag()
    {
      static util::ShPtr< command::FlagInterface> s_model_path_flag
      (
        new command::FlagStatic
        (
          "model_path",
          "change path for reading and writing models",
          command::Parameter
          (
            "model_path_name",
            "relative or absolute model path",
            command::ParameterCheckFileInSearchPath
            (
              "model",
              GetVersion().IsLicense()
              ? GetVersion().GetInstallPrefix() + "/model/"
              : BCL_MODEL_PATH,
              io::Directory::e_Dir
            ),
            ""
          )
        )
      );

      return s_model_path_flag;
    }

    //! @brief enum, that adds model params flag to default app flags
    static const util::ShPtr< command::FlagInterface> e_ModelPathFlag
    (
      command::GetAppDefaultFlags().AddDefaultFlag( Model::GetModelPathFlag(), command::e_Model)
    );

    //! @brief identifier for the name space
    //! @return the name of the namespace
    const std::string &GetNamespaceIdentifier()
    {
      static const std::string *s_namespace_name( new std::string( util::ExtractNamespaceIdentifier( __PRETTY_FUNCTION__)));
      return *s_namespace_name;
    }

    //! @brief given a FILENAME, the model path is prepended to the filename
    //! @param FILENAME the filename to a model that is used in one of the scores
    //! @return string with modelpath/FILENAME
    std::string Model::AddModelPath( const std::string &FILENAME)
    {
      const std::string resolved_filename
      (
        util::GetRuntimeEnvironment().ResolveFileName
        (
          GetModelPathFlag()->GetFirstParameter()->GetValue() + FILENAME
        )
      );

      BCL_Assert
      (
        !resolved_filename.empty(),
        "unable to resolve filename: " + GetModelPathFlag()->GetFirstParameter()->GetValue() + FILENAME
      );

      return resolved_filename;
    }

  } // namespace model
} // namespace bcl
