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
#include "command/bcl_command_parameter_check_interface.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_validation_result.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace command
  {

    //! @brief process the parameter, returning ValidationResult, and cleaning the parameter if it was valid
    //! @param PARAMETER the parameter to process
    //! @param PARAMETER_NAME the name of the parameter being processed
    //! @param ERROR_STREAM the stream to which errors are written
    //! @return returns ValidationResult
    io::ValidationResult ParameterCheckInterface::Check
    (
      const std::string &PARAMETER,
      const std::string &PARAMETER_NAME,
      std::ostream &ERROR_STREAM
    ) const
    {
      // by default, check whether the parameter is equal to the global help string, if so, write the parameter's help
      // and return true
      if( PARAMETER == io::ValidationResult::GetHelpString())
      {
        ERROR_STREAM << "Requirements for <" << PARAMETER_NAME << ">:\n";
        WriteHelp( ERROR_STREAM);
        return io::e_Help;
      }

      // check whether the parameter was allowed
      return io::ValidationResult( IsAllowedParameter( PARAMETER, PARAMETER_NAME, ERROR_STREAM));
    }

  } // namespace command
} // namespace bcl
