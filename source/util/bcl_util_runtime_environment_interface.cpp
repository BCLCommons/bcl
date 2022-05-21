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
#include "util/bcl_util_runtime_environment_interface.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_runtime_environments.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace util
  {

  ////////////////
  // operations //
  ////////////////

    //! @brief Returns the operating system specific path separator
    //! @return the operating-system specific file path separator
    const std::string &RuntimeEnvironmentInterface::GetPathSeperator() const
    {
      static const std::string s_path_seperator( "/");

      // return
      return s_path_seperator;
    }

    //! @brief report progress
    //! @param FRACTION_PROGRESS progress between [0,1]
    //! @return true if successful reported
    bool RuntimeEnvironmentInterface::ReportProgress( const double FRACTION_PROGRESS) const
    {
      BCL_MessageStd( Format().W( 3).FFP( 0)( 100 * FRACTION_PROGRESS) + "% finished");
      return true;
    }

    //! @brief add name to the process number
    //! @param NAME name to be added to process number
    //! @return NAME + processnumber as string
    std::string RuntimeEnvironmentInterface::AddProcessNumberToName( const std::string &NAME) const
    {
      return std::string( NAME + Format()( GetProcessNumber()));
    }

    //! @brief returns active RuntimeEnvironment
    //! @return active RuntimeEnvironment
    const RuntimeEnvironmentInterface &GetRuntimeEnvironment()
    {
      return GetRuntimeEnvironments().GetCurrentEnvironment();
    }

  } // namespace util
} // namespace bcl
