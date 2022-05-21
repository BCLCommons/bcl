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
#include "restraint/bcl_restraint_handler_interface.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_flag_dynamic.h"
#include "io/bcl_io_directory_entry.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {
    //! @brief get the restraints filename
    //! @return filename for this restraints; uses -restraint_prefix flag and GetFilenameExtension()
    std::string HandlerInterface::GetFilename() const
    {
      return GetFlagRestraintsFilePrefix()->GetFirstParameter()->GetValue() + GetFilenameExtension();
    }

    //! @brief test for the existence of the restraints file
    bool HandlerInterface::Exists() const
    {
      return io::DirectoryEntry( GetFilename()).DoesExist();
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer HandlerInterface::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription
      (
        "Allows use of " + GetAlias() + " restraints. The restraints are read in from X."
        + m_DefaultExtension + ", where X is the path/prefix given to -restraint_prefix "
      );
      if( m_DefaultExtension.empty())
      {
        serializer.AddInitializer( "extension", "File extension for this restraint type. The restraints are read in from {path/prefix given to -restraint_prefix}.{extension}", io::Serialization::GetAgent( &m_Extension));
      }
      else
      {
        serializer.AddInitializer
        (
          "extension",
          "File extension, defaults to " + m_DefaultExtension
          + ". The restraints are read in from {path/prefix given to -restraint_prefix}.{extension}",
          io::Serialization::GetAgent( &m_Extension),
          ""
        );
      }
      return serializer;
    }

  } // namespace restraint
} // namespace bcl
