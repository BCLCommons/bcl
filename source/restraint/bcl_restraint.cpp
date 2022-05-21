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
#include "restraint/bcl_restraint.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_flag_dynamic.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter.h"
#include "command/bcl_command_parameter_check_serializable.h"
#include "restraint/bcl_restraint_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {

    //! @brief return command line flag for specifying desired restraint types to use
    //! @return command line flag for specifying desired restraint types to use
    const util::ShPtr< command::FlagInterface> &GetFlagRestraintsTypes()
    {
      // initialize static flag
      static util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagDynamic
        (
          "restraint_types",
          "\tone or more restraint types to be used",
          command::Parameter
          (
            "restraint_type",
            "any restraint type from the list",
            command::ParameterCheckSerializable( Type())
          ),
          0,
          util::Enumerated< Interface>::GetSize()
        )
      );

      return s_flag;
    }

    //! @brief return command line flag for specifying the prefix prepended to each restraint's postfix
    //! @return command line flag for specifying the prefix prepended to each restraint's postfix
    util::ShPtr< command::FlagInterface> &GetFlagRestraintsFilePrefix()
    {
      // initialize static flag
      static util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic
        (
          "restraint_prefix", "\tthe prefix prepended to each restraint's postfix",
          command::Parameter( "prefix", "\tthe prefix prepended to each restraint's postfix", "")
        )
      );
      // end
      return s_flag;
    }

    //! @brief identifier for the name space
    //! @return the name of the namespace
    const std::string &GetNamespaceIdentifier()
    {
      static const std::string *s_namespace_name( new std::string( util::ExtractNamespaceIdentifier( __PRETTY_FUNCTION__)));
      return *s_namespace_name;
    }

  } // namespace restraint
} // namespace bcl
