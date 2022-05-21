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

// initialize the initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "io/bcl_io_serialization.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_parameter_check_file_existence.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace io
  {

    //! @brief handler for input files; checks that the file exists before setting
    //! @param VALUE a pointer to the data member to obtain a handler for
    //! @return a serialization handler for the given object
    util::OwnPtr< SerializationBase< std::string> >
      Serialization::GetAgentInputFilename( const std::string *VALUE)
    {
      return
        util::OwnPtr< SerializationBase< std::string> >
        (
          new SerializationWithCheck< std::string>( command::ParameterCheckFileExistence(), VALUE)
        );
    }

  } // namespace io
} // namespace bcl
