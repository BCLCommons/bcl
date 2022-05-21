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

#ifndef BCL_COMMAND_PARAMETER_CHECK_INTERFACE_H_
#define BCL_COMMAND_PARAMETER_CHECK_INTERFACE_H_

// include the namespace header
#include "bcl_command.h"

// include other forward headers - sorted alphabetically
#include "io/bcl_io.fwd.hh"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_object_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace command
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ParameterCheckInterface
    //! @brief This class is an interface for classes that check command line parameters
    //!
    //! @remarks example unnecessary
    //! @author karakam, woetzen, heinzes1, mendenjl
    //! @date 09/11/2009
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ParameterCheckInterface :
      public util::ObjectInterface
    {

    public:

    ////////////////
    // operations //
    ////////////////

      //! @brief returns if PARAMETER is an allowed parameter. overwritten by derived classes
      //! @param PARAMETER the parameter to check
      //! @param PARAMETER_NAME the name of the parameter being checked
      //! @param ERROR_STREAM the stream to which errors are written
      //! @return returns true if parameter is allowed, false otherwise
      virtual bool IsAllowedParameter
      (
        const std::string &PARAMETER,
        const std::string &PARAMETER_NAME,
        std::ostream &ERROR_STREAM
      ) const = 0;

      //! @brief process the parameter, returning ValidationResult, and cleaning the parameter if it was valid
      //! @param PARAMETER the parameter to process
      //! @param PARAMETER_NAME the name of the parameter being processed
      //! @param ERROR_STREAM the stream to which errors are written
      //! @return returns ValidationResult
      virtual io::ValidationResult Check
      (
        const std::string &PARAMETER,
        const std::string &PARAMETER_NAME,
        std::ostream &ERROR_STREAM
      ) const;

      //! @brief clean the parameter value; that is, make it usable by code that is independent of the parameter check
      //! Some checks may include or-type functionality, such as checking for a file among a list of paths
      //! The clean function updates the parameter value with the correct path
      //! @param PARAMETER the parameter to clean
      virtual void Clean( std::string &PARAMETER) const
      {
      }

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief writes the help for the commandline
      //! @param OSTREAM output stream to be written to
      //! @param INDENT amount to indent each line after the first
      //! @return ostream which was written to
      virtual std::ostream &WriteHelp( std::ostream &OSTREAM, const size_t INDENT = 0) const = 0;

    }; // class ParameterCheckInterface

  } // namespace command
} // namespace bcl

#endif //BCL_COMMAND_PARAMETER_CHECK_INTERFACE_H_
