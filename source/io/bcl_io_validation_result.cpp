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
#include "io/bcl_io_validation_result.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_command_state.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace io
  {
  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! constructor from bool; appropriate when the parameter is never help
    ValidationResult::ValidationResult( const bool &VALID) :
      m_Type( VALID ? e_Allowed : e_Invalid)
    {
    }

    //! constructor from Type
    ValidationResult::ValidationResult( const ValidationResultType &TYPE) :
      m_Type( TYPE)
    {
      if( TYPE == e_Help && command::CommandState::GetInMainCommandLineParsing())
      {
        command::CommandState::GetWasHelpRequested() = true;
        // only in exceptional cases is help requested but not immediately provided by the flag itself.
        // in these cases, it is expected for the application to handle the help; so the verbose command line settings
        // are unnecessary and distracting to the user then
        command::CommandState::GetWasHelpGiven() = true;
      }
    }

  //////////////////////
  // input and output //
  //////////////////////

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief get the string the user can use to indicate help is desired for a particular parameter
    //! @return help
    const std::string &ValidationResult::GetHelpString()
    {
      static const std::string s_help( "help");
      return s_help;
    }
  } // namespace io
} // namespace bcl
